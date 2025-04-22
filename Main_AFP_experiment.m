%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytic Fourier Ptychotomography (AFP) experimental data reconstruction code 
% Written by Zhenyu Dong, Haowen Zhou and Ruizhi Cao
% @ Caltech Biophotonics Laboratory | http://biophot.caltech.edu/
% The source code is licensed under GPL-3. 
% Version: March, 20th, 2025
%
% Reference:
% Project Page:
% GitHub Repo:
% 
% Organization of the code (Steps):
% 1. Set reconstruction parameters and options
% 2. Load experimental data and preprocessing
% 3. Set illumination angles / NA
% 4. AFP reconstruction
    % 4-1 Analytic NA-matching spectrum recovery using K-K relation
    % 4-2 Analytic aberration correction 
    % 4-3 Stitch NA-matching sub-spectra
    % 4-4 Analytic darkfield reconstruction (Optional) 
% 5. Generate absorption images (Optional)
% 6. Save Results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How to run the code for your own data?
% 1. Prepare the data in the same format as the example data (.mat file)
  % Data File (.mat) variables:
    % imStack: 3D data stack [NumOfPixels along x,NumOfPixels along y,numImages]
        % format: "uint8", "uint16", "double" or "single"
    % ExposureTime: exposure time for each image (normalization purpose)
    % NA_Illu: illumination angle list (N-by-2).
        % Note 1: illumination_NA = sqrt(NA_Illu(idx,1)^2 + NA_Illu(idx,1)^2)
        % Note 2: NA_Illu is dimensionless and the values belongs to [0,1]
% 2. Set the sample name and data path in Step 1
% 3. Set the parameters in setSystemParameters.m (design you own options)
% 4. Add calculation of "kIllu" in Step 3 in your sample name option
    % Code example:
    % kIllu = NA_Illu / NA * maxCTF; % unit: 1/um
    % kIllu(1:nNAmatching,:) = kIllu(1:nNAmatching,:)*0.98; % for NA matching
% 5. Run the code

% Please refere to "% May need to be modified if you have your own data"
% in the code for modifying the code for your own data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; close all; clc;
%% Check dependencies and add paths
% Package requirements:
% 1. Image Processing Toolbox
% 2. Symbolic Math Toolbox
% 3. Parallel Computing Toolbox (Optional, for GPU computing)
if ~any(any(contains(struct2cell(ver), 'Image Processing Toolbox')))
    error('Image processing toolbox is required.');
end
if ~any(any(contains(struct2cell(ver), 'Symbolic Math Toolbox')))
    error('Symbolic math toolbox is required.');
end
if ~any(any(contains(struct2cell(ver), 'Parallel Computing Toolbox'))) && (gpuDeviceCount("available") ~= 0)
    error('Parallel Computing toolbox is required.');
end

addpath(genpath('subfunctionAFP/'));

%% ------------------------------------------------------------------------
% Step 1: Set reconstruction parameters and options
% -------------------------------------------------------------------------
sampleName          = "bead";   % sample name: "bead", "embryo", "root" or "pathology"
approxMethod        = "Born";   % weak scattering approximation method: "Born" or "Rytov"
useAbeCorrection    = true;     % whether to correct for aberration in the reconstruction
saveResults         = true;     % whether to save the reconstruction results
useGPU              = true;     % whether to use GPU for reconstruction 
                                % (automatic set to false if no GPU is available)
useDarkfield        = true;     % whether to use darkfield reconstruction
enableROI           = false;    % whether to select a specific ROI for reconstruction 
                                % (set to false for all the example data)

if gpuDeviceCount("available") == 0
    useGPU = false;
end

% -------- Add data file here !!! ---------------------------------------
dataRootPath        = './data';                % data root path
dataFolder          = 'Bead';                  % data folder name
dataName            = 'Data_Bead_Upright';     % data name

dataPath            = fullfile(dataRootPath,dataFolder,[dataName,'.mat']); 

% Check if the data path is correct
if exist(dataPath, 'file')~=2
    error('Wrong data path!');
end

% Set reconstruction parameters
nNAmatching = 48;
% -------------------------------------------------------------------------
% May need to be modified if you have your own data
p = setSystemParameters(sampleName);
% -------------------------------------------------------------------------
names = fieldnames(p);
for i = 1:numel(names)
    assignin('caller', names{i}, p.(names{i}));
end

%% ------------------------------------------------------------------------ 
% Step 2: Load experimental data and preprocessing
% -------------------------------------------------------------------------
% Load experimental data
load(dataPath);

if enableROI
    ROISize = input(['Enter your ROI size, maximum value is ',num2str(min(size(imStack,1),size(imStack,2))),' : ']);
    I_sum = sum(double(imStack(:,:,1:nNAmatching)),3);
    fig = figure('Name','Incoherent image'); imshow(I_sum,[]);
    h = imrect(gca,[10 10 ROISize ROISize]);
    position = wait(h);
    position = round(position);
    close(fig);
    imStack = imStack(position(2):position(2)+position(4)-1,position(1):position(1)+position(3)-1,:);
end

imStack = double(imStack);

xsize = size(imStack,1);
ysize = size(imStack,2);  % prefer to be the same as xsize

% (1) Exposure time normalization
% -------------------------------------------------------------------------
% May need to be modified if you have your own data
imStack_temp = imStack;
imStack_temp(:,:,1:nNAmatching) = imStack(:,:,1:nNAmatching) .* ExposureTime(end) / ExposureTime(1);
for r = 1:length(ExposureTime) -1
    startFrame = nNAmatching + 1 + (sum(nDFperR(1:r)) - nDFperR(r));
    endFrame = nNAmatching + sum(nDFperR(1:r));
    imStack_temp(:,:,startFrame:endFrame) = imStack(:,:,startFrame:endFrame) .* ExposureTime(end) / ExposureTime(r+1);
end
% -------------------------------------------------------------------------

% (2) Normalize NA-matching measurements to have mean values of 1 
meanval = mean2(imStack_temp(:,:,1:nNAmatching));
imStack_temp = imStack_temp ./ meanval;
imStack = imStack_temp;
clear imStack_temp

% Generate Coherent Transfer Function (CTF)
res_xy = 1/(xsize*pixelSizeXY);           % pixel size in x-y frequency domain
res_z = 1/(zsize*pixelSizeZ);             % pixel size in z frequency domain
maxCTF = round(NA/(lambda*10^-3)/res_xy); % calculate the maximal spatial frequency in FT domain
[Y,X] = meshgrid(1:ysize,1:xsize);
xc = floor(xsize/2+1);
yc = floor(ysize/2+1);
R_k4img = abs((X-xc) + 1i*(Y-yc));
CTF = R_k4img<= maxCTF;

%% ------------------------------------------------------------------------
% Step 3: Set illumination angles / NA
% -------------------------------------------------------------------------
% kIllu is a N by 2 matrix, containing the (kx,ky) position in pixels for all N angles
% N (only for illustration here) is the number of illumination angles (NA-matching + darkfield)

% -------------------------------------------------------------------------
% May need to be modified if you have your own data
if sampleName == "bead"
    % na_calib is the illumination NA obtained from self-calibration
    kIllu = fliplr(na_calib/(lambda*10^-3)/res_xy);  % unit: pixel

elseif sampleName == "embryo"
    % NA_Illu_scan is the designed illumination NA set in galvo mirror
    % scanning code, therefore it remains a circle shape, and we fit an
    % elliptical function for the voltage values.
    NA_Illu_scan(1:nNAmatching,:) = NA_Illu_scan(1:nNAmatching,:)*0.98;
    kIllu = fliplr(NA_Illu_scan/NA*maxCTF);          % unit: pixel

elseif sampleName == "root"
    % Eight rings in total, since root is thick and requires more overlaps
    % for the darkfield reconstruction
    phi = pi/56;
    myAngles = linspace(0,2*pi,nNAmatching+1)+phi;
    kIllu(:,1) = cos(myAngles(1:end-1))*maxCTF*0.98;
    kIllu(:,2) = -sin(myAngles(1:end-1))*maxCTF*0.98;
    
    NA_list = [0.4018, 0.4458, 0.4735, 0.5057, 0.5358, 0.5640, 0.5774];
    NmaxCTF =  round(NA_list./(lambda*10^-3)/res_xy);
    for idring = 2:floor(size(imStack,3)/nNAmatching)
        kIllu((idring-1)*nNAmatching+1:idring*nNAmatching,1) = cos(myAngles(1:end-1))*NmaxCTF(idring);
        kIllu((idring-1)*nNAmatching+1:idring*nNAmatching,2) = -sin(myAngles(1:end-1))*NmaxCTF(idring);
    end

elseif sampleName == "pathology"
    % Since a single UV LED rotates with a rotation stage, the angles lie
    % on a ring shape
    phi = 6/180*pi;
    myAngles = (myAngles + phi);
    kIllu(:,1) = sin(myAngles(1:end-1))*maxCTF*0.97;
    kIllu(:,2) = cos(myAngles(1:end-1))*maxCTF*0.97;

else
    disp('Illumination angles for this sample is not defined');
    disp('sampleName should be : "bead", "embryo", "root" or "pathology"')
    disp('Or add your own experimental parameters.');
    disp(' ')
end
% -------------------------------------------------------------------------


p.kIllu = kIllu;

figure('Name','Illumination angles');
plot(kIllu(:,1)/maxCTF*NA,kIllu(:,2)/maxCTF*NA,'o:');
xlabel('NA_x'); ylabel('NA_y');
axis equal; axis tight;
title('Illumination angle trajectory');

%% ------------------------------------------------------------------------
% Step 4: AFP reconstruction 
% Step 4.1: Analytic NA-matching spectrum recovery using K-K relation
% -------------------------------------------------------------------------
[recFTframe,mask2use] = recFieldKK(imStack(:,:,1:nNAmatching),kIllu(1:nNAmatching,:),...
                                   'CTF',CTF,'pad',4,'norm',true,'wiener',true); 
                        % Outputs:
                        % recFTframe: Sub-spectra of reconstructed total field (unscattered+scattered) under NA-matching measurements
                        % mask2use: A binary mask for CTF

% Check illumination angle order, should have the same position with DC in recFTframe
% scatter(kIllu(1:8,2),-kIllu(1:8,1));xlim([-maxCTF,maxCTF]);ylim([-maxCTF,maxCTF]);axis square;

%% ------------------------------------------------------------------------ 
% Step 4: AFP reconstruction
% Step 4-2: Analytic aberration correction
% -------------------------------------------------------------------------
k_carrier = n_media/(lambda*10^-3); % unit: 1/um. wavenumber of the incident light
kzMap = (sqrt(k_carrier^2 - (R_k4img*res_xy).^2)).*(abs(CTF)>10^-3); 
% unit: 1/um. z-axis wavenumber for each point in the Fourier domain, it has the same shape as the Ewald's sphere. 

if useAbeCorrection
    % CTF_abe is the retrieved aberration function
    % zernikeCoeff is zernike coeffient starting from piston term
    [CTF_abe,zernikeCoeff] = findAbeFromOverlap3D(recFTframe,kIllu(1:nNAmatching,:),pixelSizeZ,CTF,kzMap,...
                                                  'weighted',true,'zernikemode',3:25,'reg',0.1, ...
                                                  'thickness',samThick,'similarity tolerence',sqrt(2)/2, ...
                                                  'vis_samplingmap',false); 
    disp('AFP aberration correction: done!');
    
    % Visualize aberration 
    cmapAberration = cNeoAlbedo;  % load colormap for aberration
    figure('Name','Retrieved aberration ');
    subplot(121); imagesc(angle(CTF_abe).*(abs(CTF_abe)>10^-2)); title('reconstructed aberration');
    clim([-pi pi]); colormap(gca,cmapAberration); axis image; axis off; colorbar;

    subplot(122); bar(zernikeCoeff);legend('AFP'); title('zernike coefficients');
    set(gcf, 'Position', [400, 400, 1000, 300]);

end

%% ------------------------------------------------------------------------
% Step 4: AFP reconstruction 
% Step 4-3: Stitching NA-Matching sub-spectra
% -------------------------------------------------------------------------
if useDarkfield
    highresPad = 2;                 % prepare for spectrum extention
else
    highresPad = 1;
end
xsizeR = xsize*highresPad;
ysizeR = ysize*highresPad;
zsizeR = zsize*highresPad;          % keep the same dimensions for isotropic padding
recVolFT = zeros(xsizeR,ysizeR,zsizeR); 
recScatteringPotentialFT = zeros(xsizeR,ysizeR,zsizeR);
weightMtx = zeros(size(recVolFT));
[YR,XR] = meshgrid(1:ysizeR,1:xsizeR);
xcR = floor(xsizeR/2+1);
ycR = floor(ysizeR/2+1);
zcR = floor(zsizeR/2+1);

% Calculate Ewald's Sphere 
kzNormMap = (sqrt((n_media/(lambda*10^-3))^2 - (R_k4img*res_xy).^2) - n_media/(lambda*10^-3))/res_z; % in pixels
kzNormMap = kzNormMap.*(abs(CTF)>10^-3);

% Stitching
indexingCTF = abs(CTF)>10^-3;
for idx = 1:nNAmatching
    k2use = round(kIllu(idx,:));
    kzOffset = (n_media/(lambda*10^-3) - sqrt((n_media/(lambda*10^-3))^2 - (kIllu(idx,1)*res_xy)^2- ...
               (kIllu(idx,2)*res_xy)^2))/res_z; % in pixels
    kz2use = round(kzNormMap + kzOffset) + zcR;
    linearIdx = sub2ind([xsizeR,ysizeR,zsizeR],X(indexingCTF)-k2use(1) - xc + xcR,...
                                               Y(indexingCTF)-k2use(2) - yc + ycR,...
                                               kz2use(indexingCTF));

    weightMtx(linearIdx) = weightMtx(linearIdx) + mask2use(indexingCTF);
    
    % Correct aberration for each sub-spectrum before stitching
    if useAbeCorrection
        offsetPhase = angle(CTF_abe(xc+round(kIllu(idx,1)),yc+round(kIllu(idx,2))));
        temp = recFTframe(:,:,idx).*conj(CTF_abe)*exp(1i*offsetPhase)./(abs(CTF_abe) + 10^-3);
    else
        temp = recFTframe(:,:,idx).*mask2use./(mask2use + 10^-3);
    end
    
    recVolFT(linearIdx) = recVolFT(linearIdx) + temp(indexingCTF);
    
    % Use Rytov/Born approximation to convert total field from K-K 
    % reconstruction into first order scattered field
    if approxMethod == "Rytov" 
        powerRatio = sqrt(sum(sum(temp.*conj(temp))))/(xsize*ysize);
        temp = ifft2(ifftshift(circshift(temp,-round(kIllu(idx,:)))))/powerRatio;
        phs = phase_unwrapCG(angle(temp));
        amp = abs(temp);
        temp = circshift(fftshift(fft2((log(amp)+1j*phs))),round(kIllu(idx,:)));
    elseif approxMethod == "Born"
        temp(xc + round(kIllu(idx,1)),yc + round(kIllu(idx,2))) = ...
            temp(xc + round(kIllu(idx,1)),yc + round(kIllu(idx,2))) - sqrt(sum(sum(temp.*conj(temp))));
    end

    % Convert first order scattered field to scattering potential using
    % Fourier diffraction theorem
    kz = (n_media/(lambda*10^-3))^2 - (R_k4img*res_xy).^2; kz(kz<0) = 0; kz = sqrt(kz);
    temp = temp*-2j.*kz;

    % Stitch the scattering potential Fourier spectrum
    recScatteringPotentialFT(linearIdx) = recScatteringPotentialFT(linearIdx) + temp(indexingCTF)*(2*pi*highresPad^3/pixelSizeZ);

end

% Convert scattering potential back to refractive index (RI_3D)
recVolFT = recVolFT./(weightMtx + 10^-5);
fieldRec = fftshift(ifftn(ifftshift(recVolFT)),3);
maskRecons = (weightMtx>10^-3);
recScatteringPotentialFT = recScatteringPotentialFT./(weightMtx + 10^-5);
ScatteringPotentialRec = fftshift(ifftn(ifftshift(recScatteringPotentialFT)),3);

% RI_3D is AFP reconstruction result using NA-matching measurements
[RI_3D,~] = convertScatteringPotentialToRI(ScatteringPotentialRec,lambda*1e-3,n_media);
disp('AFP NA-matching reconstruction: done!');

%% ------------------------------------------------------------------------
% Step 4: AFP reconstruction
% Step 4-4: Analytic darkfield reconstruction
% -------------------------------------------------------------------------
if useDarkfield
    % Note: (1) 'reg' is used to consider the misalignment of illumination
    % angles, which usually occurs in the experiment setups. Therefore, a
    % larger 'reg' can give a better result that is free from grainy
    % artifact. Generally, set 'reg' as 2 in experiments, no need to tune.
    % 
    % Note: (2) mask_fullfield contains: NA-matching + entire darkfield, while
    % maskRecons_fullfield contains: NA-matching + unknown part of darkfield
    [recVolFT,maskRecons_fullfield,recScatteringPotentialFT_fullfield,mask_fullfield] = recFieldFromKnown3D( ...
                                                 imStack(:,:,nNAmatching+1:end),kIllu(nNAmatching+1:end,:),...
                                                 recVolFT,maskRecons,recScatteringPotentialFT,CTF_abe,k_carrier, ...
                                                 kzNormMap,indexingCTF,zcR,res_xy,res_z,kz,pixelSizeZ,...
                                                 'drift',true,'reg',2,'approx',approxMethod,'gpu',useGPU,...
                                                 'intensity correction', true, 'thres', 0.30,...
                                                 'use data intensity',true,'thickness',samThick, ...
                                                 'similarity tolerence',0.3,'timer on','darkfield');

    % fieldRecFull = fftshift(ifftn(ifftshift(recVolFT)),3); % This is total field

    % Scattering potential (NA-matching + entire darkfield)
    ScatteringPotentialRec_fullfield = fftshift(ifftn(ifftshift(recScatteringPotentialFT_fullfield)),3);
    % Convert back to refractive index
    [RI_3D_fullfield,~] = convertScatteringPotentialToRI(ScatteringPotentialRec_fullfield,lambda*1e-3, n_media);
end

%% ------------------------------------------------------------------------ 
% Step 5: Generate absorption images (Optional) 
% -------------------------------------------------------------------------
if useAbsorption
    % Real and imaginary part of the refractive index
    [RI_3D,Absorb_3D] = convertScatteringPotentialToRI(ScatteringPotentialRec,lambda*1e-3,n_media);
    
    % Complex-valued refractive index
    RI_recons = RI_3D+1j*Absorb_3D;
    padSize = 0; 

    % Synthesize absorption volume 'imStack_ICNAMatch' by simulating 
    % incoherent illumination
    imStack_ICNAMatch = zeros(size(RI_recons));
    for idx_illu = 1:nNAmatching
        if mod(idx_illu,10)==0
            disp(['angle: ',num2str(idx_illu)]);
        end
        % Since we can separte out imaginary part of the refractive index, 
        % we emphasize the absorption effect by assuming that real part is 
        % the same as the background medium.
        imStack_temp = imagingMultiSliceStack(n_media+1j*imag(RI_recons),kIllu(idx_illu,:),CTF,lambda, ...
                                              pixelSizeXY,pixelSizeZ,'n',n_media,'pad',padSize,'gpu',useGPU);
        imStack_ICNAMatch = imStack_ICNAMatch + imStack_temp;
    end

    figure('Name','AFP reconstruction');
    subplot(121);imshow(RI_3D(:,:,floor(zsize/2+1)),[]);clim([1.4975,1.502]);colorbar;title('AFP RI real part');
    subplot(122);imshow(mat2gray(imStack_ICNAMatch(:,:,floor(zsize/2+1)+5)),[]);colorbar;title('AFP absorption');
    
    % Put the complex RI back to the forward model and compare the predicted 
    % measurements with the real measurements on the camera
    imStack_predict = imagingMultiSlice(RI_recons,nNAmatching,kIllu,CTF,lambda,pixelSizeXY,pixelSizeZ, ...
                                        'n',n_media,'pad',padSize,'gpu',useGPU,'use_darkfield',false);
    figure('Name','AFP Absorption reconstruction validation');
    subplot(121);imshow(sum(imStack,3),[]);colorbar;title('GT Measurement (Incoherent)');
    subplot(122);imshow(sum(imStack_predict,3),[]);colorbar;title('AFP Prediction (Incoherent)');
end

%% ------------------------------------------------------------------------
% Step 6: Save Results
% -------------------------------------------------------------------------
if saveResults
    % save reconstruction parameters
    p.sampleName = sampleName;
    p.approxMethod = approxMethod;   
    p.useAbeCorrection = useAbeCorrection; 
    p.useDarkfield = useDarkfield;
    p.highresPad = highresPad;

    if useAbeCorrection
        CTF_abe = CTF_abe.*(abs(CTF_abe)>10^-2);
    else
        CTF_abe = CTF;
        zernikeCoeff = [];
    end

    saveDir = 'Result_experiment';
    if ~exist(saveDir,'dir')
       mkdir(saveDir);
    end
   
    if useDarkfield
        file_name = fullfile(saveDir,[dataName(find(dataName == '_', 1) + 1:end),'_AFP_experiment_',char(approxMethod),'_','Darkfield','.mat']);
        save(file_name,'RI_3D_fullfield','RI_3D','CTF_abe','zernikeCoeff','p','-v7.3');
    elseif useAbsorption
        file_name = fullfile(saveDir,[dataName(find(dataName == '_', 1) + 1:end),'_AFP_experiment_',char(approxMethod),'_','NAmatching','.mat']);
        save(file_name,'RI_3D','imStack_ICNAMatch','CTF_abe','zernikeCoeff','p','-v7.3');
    else
        file_name = fullfile(saveDir,[dataName(find(dataName == '_', 1) + 1:end),'_AFP_experiment_',char(approxMethod),'_','NAmatching','.mat']);
        save(file_name,'RI_3D','CTF_abe','zernikeCoeff','p','-v7.3');
    end
end

