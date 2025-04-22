%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytic Fourier Ptychotomography (AFP) simulation code 
% Written by Ruizhi Cao, Zhenyu Dong, and Haowen Zhou
% @ Caltech Biophotonics Laboratory | http://biophot.caltech.edu/
% The source code is licensed under GPL-3. 
% Version: March, 15th, 2025
%
% Reference:
% Project Page:
% GitHub Repo:
% 
%
% Organization of the code:
% 1. Set imaging system and reconstruction parameters
    % 1-1 Simulation Options
    % 1-2 Imaging System Parameters
    % 1-3 Sample Configuration
    % 1-4 Generate Coherent Transfer Function (CTF) with aberration
    % 1-5 Design Illumination angles (NA-matching & darkfield)
% 2. Simulate raw intensity measurements
% 3. AFP reconstruction
    % 3-1 Analytic NA-matching spectrum recovery using K-K relation
    % 3-2 Analytic aberration correction (sub-spectrum stitching included)
    % 3-3 Analytic darkfield reconstruction (Optional)
% 4. Visualization & comparison with ground truth
% 5. Save results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;close all;clc;

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
addpath(genpath('sampleFunctions/'));

%% ------------------------------------------------------------------------
% Step 1: Set imaging system and reconstruction parameters
% Setp 1-1:  Simulation Options
% -------------------------------------------------------------------------
sampleName = "bead";      % sample type: "bead" or "lymphn"
approxMethod = "Born";    % weak scattering approximation method: "Born" or "Rytov"
useAbe = true;            % whether to introduce aberration in the forward model
useAbeCorrection = true;  % whether to correct for aberration in the reconstruction
usePadding = true;        % whether to pad the object in multi-slice simulation of raw measurements
saveResults = true;       % whether to save the reconstruction results
useGPU = true;            % whether to use GPU for reconstruction (automatic set to false if no GPU is available)
useDarkfield = false;     % whether to use darkfield reconstruction

if gpuDeviceCount("available") == 0
    useGPU = false;
end

%% ------------------------------------------------------------------------
% Step 1: Set imaging system and reconstruction parameters
% Setp 1-2:  Imaging System Parameters
% -------------------------------------------------------------------------
mag         = 20;         % magnification of the objective lens
ps          = 6.5/mag;    % unit: um. Pixel size in object space
pixelSizeXY = ps;         % unit: um. Effective pixel size (lateral)
pixelSizeZ  = ps;         % unit: um. Allows uneven 3D grid, i.e. pixelSizeZ can be different from pixelSizeXY
NA          = 0.4;        % numerical aperture (NA) of the imaging system
lambda      = 532;        % unit: nm. wavelength of the illumination light


%% ------------------------------------------------------------------------
% Step 1: Set imaging system and reconstruction parameters
% Setp 1-3:  Sample Configuration
% -------------------------------------------------------------------------
xsize = 200;              % unit: pixel
ysize = 200;              % unit: pixel. Prefer to be the same as xsize
zsize = 200;              % unit: pixel. zsize can be different from xsize and ysize

if sampleName == "lymphn"
    % Lymph node sample
    r2useLib    = [1.11, 1.20, 1.28, 1.37, 1.45, 1.5]; 
                              % darkfield illumination angles, illumination NA = r2useLib.*NA
    n_media     = 1.33068;    % refractive index of the background environment
    samThick    = 40;         % sample thickness in pixel before highresPad (compared to zsize)
    obj         = lymph_node([xsize,ysize,zsize],n_media,samThick);

elseif sampleName == "bead"
    % Bead sample      
    r2useLib     = [1.15, 1.28, 1.4, 1.5]; 
                              % darkfield illumination angles, illumination NA = r2useLib.*NA
    n_media      = 1.58;      % refractive index of the background environment
    diameterBead = 6;         % unit: micron.
    rBead        = diameterBead/2./[pixelSizeXY,pixelSizeXY,pixelSizeZ]; 
                              % unit: pixel. radius of the bead
    RI_Bead      = 1.59;      % refractive index of the bead sample
    samThick     = min(round(diameterBead/pixelSizeZ*1.5),round(50/pixelSizeZ)); 
                              % in pixels, compared to original obj size
    obj          = phantom_bead([xsize,ysize,zsize],rBead*2,n_media,RI_Bead);

else
    % Define your own sample function here for simulation
    disp('Sample name not defined, Options: bead | lymphn');
end

%% ------------------------------------------------------------------------
% Step 1: Set imaging system and reconstruction parameters
% Setp 1-4:  Generate Coherent Transfer Function (CTF) with aberration 
% -------------------------------------------------------------------------
res_xy = 1/(xsize*pixelSizeXY);           % unit: 1/um. pixel size in x-y frequency domain
res_z = 1/(zsize*pixelSizeZ);             % unit: 1/um. pixel size in z frequency domain
maxCTF = round(NA/(lambda*10^-3)/res_xy); % unit: pixel. Calculate the maximal spatial frequency in Fourier domain
[Y,X] = meshgrid(1:ysize,1:xsize);
xc = floor(xsize/2+1);
yc = floor(ysize/2+1);
R_k4img = abs((X-xc) + 1i*(Y-yc));
CTF = R_k4img<= maxCTF;                

if useAbe
    % Construct aberration
    nZernike = 25;                             % maximal number of Zernike modes
    zernikeGT = zeros(nZernike,1);
    temp2use = [0,0,-0.5,-1,0.6,1,0.7,-1];     % change zernike coeffs here, start from tilt term (2nd term)
                                               % Ref: https://en.wikipedia.org/wiki/Zernike_polynomials
    zernikeGT(1:length(temp2use)) = temp2use;  % Ground truth zernike coeffs

    temp = linspace(-1,1,2*maxCTF+1);
    [Yz,Xz] = meshgrid(temp,temp);
    [theta,r] = cart2pol(Yz,Xz);
    idx2use = (r<=1);
    zernikeTemp = zernfun2(1:nZernike,r(idx2use),theta(idx2use));
    Hz = zeros((2*maxCTF+1)^2,nZernike);       % Zernike operator that generates aberration
    for idx = 1:nZernike
        Hz(idx2use,:) = zernikeTemp;
    end
    clear zernikeTemp theta r Xz Yz;
    temp = reshape(Hz*zernikeGT,[2*maxCTF+1,2*maxCTF+1]);
    CTF_GT = double(CTF);
    bd = calBoundary([xc,yc],[2*maxCTF+1,2*maxCTF+1]);
    CTF_GT(bd(1):bd(2),bd(3):bd(4)) = CTF_GT(bd(1):bd(2),bd(3):bd(4)).*exp(1i*(temp));
    Hz_temp = Hz;
    clear temp bd Hz;
else
    CTF_GT = CTF;                              % CTF_GT is the ground truth aberration (amplitude & phase)
end

%% -------------------------------------------------------------------------
% Step 1: Set imaging system and reconstruction parameters
% Setp 1-5:  Design Illumination angles (NA-matching & darkfield)
% -------------------------------------------------------------------------
% NA-matching measurements
nNAmatching = 48;
myAngles = linspace(0,2*pi,nNAmatching+1);
kx_inCA = cos(myAngles(1:end-1))*maxCTF*0.98;  % unit: pixel
ky_inCA = sin(myAngles(1:end-1))*maxCTF*0.98;

% darkfield measurements
baseNum = 40;                                  % number of angles at radius equals 1 (NA).
nDFperR = ceil(r2useLib*baseNum);              % r2useLib is defined at the sample functions section
                                               % illumination NA = r2useLib.*NA
kIllu = zeros(nNAmatching + sum(nDFperR),2);
                % kIllu is a N by 2 matrix, containing the (kx,ky) position in pixels for all N angles
                % N (only for illustration here) is the number of illumination angles (NA-matching + darkfield)
idxSt = nNAmatching+1;
kIllu(1:nNAmatching,1) = kx_inCA;
kIllu(1:nNAmatching,2) = ky_inCA;

for idx = 1:length(r2useLib)
    r2use = r2useLib(idx)*maxCTF;
    myAngles = linspace(0,2*pi,nDFperR(idx)+1);
    kIllu(idxSt:idxSt+nDFperR(idx)-1,1) = cos(myAngles(1:end-1))*r2use;
    kIllu(idxSt:idxSt+nDFperR(idx)-1,2) = sin(myAngles(1:end-1))*r2use;
    idxSt = idxSt + nDFperR(idx);
end

figure('Name','Illumination angles');
plot(kIllu(:,1)/maxCTF*NA,kIllu(:,2)/maxCTF*NA,'o:');
xlabel('NA_x'); ylabel('NA_y');
axis equal; axis tight;
title('Illumination angle sequence');

%% -------------------------------------------------------------------------
% Step 2: Simulate raw intensity measurements
% Generate raw intensity measurements under all illumination angles
% --------------------------------------------------------------------------
if usePadding  
    % This guarantees that the unscattered (DC) light will cover the whole FOV in the entire simulation volume
    % (will result in horizontal/vertical lines in the Foureier spectrum, which are usually observed in experiments / natural images)
    pad_size = ceil(tan(asin(NA))*pixelSizeZ*zsize/pixelSizeXY)+3;   % pad first and crop at last step       
    maxCTF_pad = round(NA/(lambda*10^-3)/(1/((xsize+2*pad_size)*pixelSizeXY)));
    [YY,XX] = meshgrid(1:(ysize+2*pad_size),1:(xsize+2*pad_size));
    CTF_GT_pad = double(abs((XX-floor((xsize+2*pad_size)/2+1)) + 1i*(YY-floor((ysize+2*pad_size)/2+1)))<= maxCTF_pad);
    if useAbe
        temp = linspace(-1,1,2*maxCTF_pad+1);
        [Yz,Xz] = meshgrid(temp,temp);
        [theta,r] = cart2pol(Yz,Xz);
        idx2use = (r<=1);
        zernikeTemp = zernfun2(1:nZernike,r(idx2use),theta(idx2use));
        Hz = zeros((2*maxCTF_pad+1)^2,nZernike);                     % zernike operator that generates aberration
        for idx = 1:nZernike
            Hz(idx2use,:) = zernikeTemp;
        end
        clear zernikeTemp theta r Xz Yz idx2use;
        temp = reshape(Hz*zernikeGT,[2*maxCTF_pad+1,2*maxCTF_pad+1]);
        bd = calBoundary([floor((xsize+2*pad_size)/2+1),floor((ysize+2*pad_size)/2+1)],[2*maxCTF_pad+1,2*maxCTF_pad+1]);
        CTF_GT_pad(bd(1):bd(2),bd(3):bd(4)) = CTF_GT_pad(bd(1):bd(2),bd(3):bd(4)).*exp(1i*(temp));
        clear temp YY XX bd Hz;
    end
else
    pad_size = 0;
    CTF_GT_pad = CTF_GT;
    maxCTF_pad = maxCTF;
end

% Use multi-sclice beam propagation method to simulate raw intensity measurements on the camera
imStack = imagingMultiSlice(obj,nNAmatching,kIllu,CTF_GT_pad,lambda,pixelSizeXY,pixelSizeZ,...
                            'n',n_media,'pad',pad_size,'gpu',useGPU,'use_darkfield',useDarkfield);

disp('AFP reconstruction start!');
%% ------------------------------------------------------------------------- 
% Step 3: AFP reconstruction
% Step 3-1: Total field reconstruction using K-K relation
% -------------------------------------------------------------------------
[recFTframe,mask2use] = recFieldKK(imStack(:,:,1:nNAmatching),kIllu(1:nNAmatching,:),...
                                   'CTF',CTF,'pad',4,'norm',true,'wiener',false); 
                        % Outputs:
                        % recFTframe: Sub-spectra of reconstructed total field (unscattered+scattered) under NA-matching measurements
                        % mask2use: A binary mask for CTF

%% -------------------------------------------------------------------------
% Step 3: AFP reconstruction 
% Step 3-2: Analytic aberration correction (sub-spectrum stitching included)
% -------------------------------------------------------------------------
k_carrier = n_media/(lambda*10^-3); % unit: 1/um. wavenumber of the incident light
kzMap = (sqrt(k_carrier^2 - (R_k4img*res_xy).^2)).*(abs(CTF)>10^-3); 
                                    % unit: 1/um. z-axis wavenumber for each point in the Fourier domain, 
                                    % and kzMap has the same shape as the Ewald's sphere. 

if useAbeCorrection
    % CTF_abe is the retrieved aberration function
    % ZernikeCoeff is zernike coeffient starting from piston term (1st term)
    % Ref: http://en.wikipedia.org/wiki/Zernike_polynomials
    [CTF_abe,zernikeCoeff] = findAbeFromOverlap3D(recFTframe,kIllu(1:nNAmatching,:),pixelSizeZ,CTF,kzMap,...
                                                  'weighted',true,'zernikemode',3:25,'reg',0.1, ...
                                                  'thickness',samThick,'similarity tolerence',sqrt(2)/2, ...
                                                  'vis_samplingmap',false); 
    disp('AFP aberration correction: done!');

    % Calculate aberration differences (tilt & defocus term removed)
    zernike_diff = zernikeCoeff.'-[0;zernikeGT];
    zernike_diff([1:3,5]) = [];         % remove tilt & defocus zernike coefficient
    Hz_temp(:,[1:2,4]) = [];            % remove tilt & defocus zernike polynomials
    temp = reshape(Hz_temp*zernike_diff,[sqrt(size(Hz_temp,1)),sqrt(size(Hz_temp,1))]);
    CTF_error = double(CTF);
    bdC = calBoundary([xc,yc],[sqrt(size(Hz_temp,1)),sqrt(size(Hz_temp,1))]);
    CTF_error(bdC(1):bdC(2),bdC(3):bdC(4)) = CTF_error(bdC(1):bdC(2),bdC(3):bdC(4)).*exp(1i*temp);

    % Load colormap for aberration visualization
    cmapAberration = cNeoAlbedo;  
    figure('Name','Aberration comparison');
    subplot(221); imagesc(angle(CTF_abe).*(abs(CTF_abe)>10^-2)); title('reconstructed aberration');
    clim([-pi pi]); colormap(gca,cmapAberration); axis image; axis off; colorbar;

    subplot(222); imagesc(angle(CTF_GT)); title('GT aberration');
    clim([-pi pi]); colormap(gca,cmapAberration); axis image; axis off; colorbar;

    subplot(223); imagesc(angle(CTF_error).*(abs(CTF_abe)>10^-2)); title('Aberration error');
    colormap(gca,cmapAberration); axis image; axis off; colorbar;

    subplot(224); % zernike coefficient comparison
    zernikeCoeff = zernikeCoeff.';
    abe(:,1) = zernikeCoeff;
    zernikeGT = [0;zernikeGT];
    abe(:,2) = zernikeGT;
    bar(abe);legend('AFP','GT'); 
end

% Stitching NA-Matching sub-spectra
highresPad = 2;                     % prepare for spectrum extention
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

%% -------------------------------------------------------------------------
% Step 3: AFP reconstruction
% Step 3-3: Analytic darkfield reconstruction
% -------------------------------------------------------------------------
if useDarkfield
    % Note: (1) 'reg' is used to consider the misalignment of illumination
    % angles, so 'reg' is set as a small value (0.01) in the simulation.
    % But in the experiment, misalignment usually occurs, therefore, a
    % larger 'reg' can give a better result that is free from grainy
    % artifact. Generally, set 'reg' as 2 in experiments. 
    % 
    % Note: (2) mask_fullfield contains: NA-matching + entire darkfield, while
    % maskRecons_fullfield contains: NA-matching + unknown part of darkfield
    [recVolFT,maskRecons_fullfield,recScatteringPotentialFT_fullfield,mask_fullfield] = recFieldFromKnown3D( ...
                                                 imStack(:,:,nNAmatching+1:end),kIllu(nNAmatching+1:end,:),...
                                                 recVolFT,maskRecons,recScatteringPotentialFT,CTF_abe,k_carrier, ...
                                                 kzNormMap,indexingCTF,zcR,res_xy,res_z,kz,pixelSizeZ,...
                                                 'drift',true,'reg',0.01,'approx',approxMethod,'gpu',useGPU,...
                                                 'intensity correction', true, 'thres', 0.30,...
                                                 'use data intensity',true,'thickness',samThick, ...
                                                 'similarity tolerence',0.3,'timer on','darkfield');

    % fieldRecFull = fftshift(ifftn(ifftshift(recVolFT)),3); % This is total field

    % Scattering potential (NA-matching + entire darkfield)
    ScatteringPotentialRec_fullfield = fftshift(ifftn(ifftshift(recScatteringPotentialFT_fullfield)),3);
    
    % Convert back to refractive index
    [RI_3D_fullfield,~] = convertScatteringPotentialToRI(ScatteringPotentialRec_fullfield,lambda*1e-3, n_media);
    
    % Visualization of Darkfield-extended Fourier spectrum
    figure('Name','DF_Spectrum_XYMIP');
    DF_spectrum = log10(1+abs(recScatteringPotentialFT_fullfield));
    temp = squeeze(max(DF_spectrum,[],3));
    imshow(temp,[]);axis image;colormap('hot');
    
    disp('AFP darkfield reconstruction: done!');
end

%% -------------------------------------------------------------------------
% Step 4: Visualization & comparison with ground truth
% --------------------------------------------------------------------------

% Pad refractive index ground truth
ftGT = fftshift(fftn(obj));
bd2use = calBoundary(floor(size(RI_3D)/2+1),size(ftGT));
ftTemp = zeros(size(RI_3D));
ftTemp(bd2use(1):bd2use(2),bd2use(3):bd2use(4),bd2use(5):bd2use(6)) = ftGT;
ftGT = ftTemp;
RIGT_pad = real(ifftn(ifftshift(ftGT)))*(highresPad)^3;

% Apply spectrum mask to ground truth scattering potential and convert back to refractive index
ScatteringPotentialGT_pad = (2*pi/(lambda*1e-3))^2*((ifftn(ifftshift(ftGT))*(highresPad)^3).^2-n_media.^2);
[RIGT_pad_LP,~] = convertScatteringPotentialToRI(ifftn(ifftshift(fftshift(fftn(ScatteringPotentialGT_pad)).*maskRecons)),lambda*1e-3,n_media);

% Draw theoretical reconstruction results (GT vs GT(NA-matching) vs AFP(NA-matching))
screenSize = get(0, 'ScreenSize');
figure('Name','RI comparison (NA-matching)','Position', screenSize);
subplot(231); imshow(RIGT_pad(:,:,zcR),[]); title('GT (x-y view)'); colorbar;
subplot(232); imshow(RIGT_pad_LP(:,:,zcR),[]); title('GT NA-matching (x-y view)'); colorbar;
subplot(233); imshow(RI_3D(:,:,zcR),[]); title('AFP NA-matching (x-y view)'); colorbar;
subplot(234); imshow(squeeze(RIGT_pad(:,ycR,:)).',[]); title('GT (x-z view)'); colorbar;
subplot(235); imshow(squeeze(RIGT_pad_LP(:,ycR,:)).',[]); title('GT NA-matching (x-z view)'); colorbar;
subplot(236); imshow(squeeze(RI_3D(:,ycR,:)).',[]); title('AFP NA-matching (x-z view)');colorbar;

% Draw darkfield reconstruction results comparison
if useDarkfield   
    [RIGT_pad_LP_fullfield,~] = convertScatteringPotentialToRI(ifftn(ifftshift(fftshift(fftn( ...
                                 ScatteringPotentialGT_pad)).*mask_fullfield)),lambda*1e-3,n_media);
    
    % Visualize reconstruction results (GT vs GT(NA-matching + DF) vs AFP(NA-matching + DF))
    figure('Name','RI comparison (Fullfield)','Position', screenSize);
    subplot(231);imshow(RIGT_pad(:,:,zcR),[]);title('GT (x-y view)');colorbar;
    subplot(232);imshow(RIGT_pad_LP_fullfield(:,:,zcR),[]);title('GT NA-matching + DF (x-y view)');colorbar;
    subplot(233);imshow(RI_3D_fullfield(:,:,zcR),[]);title('AFP NA-matching + DF (x-y view)');colorbar;
    subplot(234);imshow(squeeze(RIGT_pad(:,ycR,:)).',[]);title('GT (x-z view)');colorbar;
    subplot(235);imshow(squeeze(RIGT_pad_LP_fullfield(:,ycR,:)).',[]);title('GT NA-matching + DF (x-z view)');colorbar;
    subplot(236);imshow(squeeze(RI_3D_fullfield(:,ycR,:)).',[]);title('AFP NA-matching + DF (x-z view)');colorbar;
    
    % Visualize reconstruction results (AFP(NA-matching) vs AFP(NA-matching + DF) vs GT(NA-matching + DF))
    figure('Name','RI comparison (NA-matching & Fullfield)','Position', screenSize);
    subplot(231);imshow(RI_3D(:,:,zcR),[]);title('AFP (NA-matching)');colorbar;
    subplot(232);imshow(RI_3D_fullfield(:,:,zcR),[]);title('AFP (NA-matching + DF)');colorbar;
    subplot(233);imshow(RIGT_pad_LP_fullfield(:,:,zcR),[]);title('GT (NA-matching + DF)');colorbar;
    subplot(234);imshow(squeeze(RI_3D(:,ycR,:)).',[]);title('AFP (NA-matching)');colorbar;
    subplot(235);imshow(squeeze(RI_3D_fullfield(:,ycR,:)).',[]);title('AFP (NA-matching + DF)');colorbar;
    subplot(236);imshow(squeeze(RIGT_pad_LP_fullfield(:,ycR,:)).',[]);title('GT (NA-matching + DF)');colorbar;

end

if useDarkfield
    RI_3D_temp = RI_3D_fullfield;
    title_temp = 'NA-matching + DF';
else
    RI_3D_temp = RI_3D;
    title_temp = 'NA-matching';
end

% Draw maximum intensity projection (MIP) rendering
fig = figure('Name','Volume rendering (MIP)'); 
% Check MATLAB version and use the appropriate volshow syntax
warnState = warning;             % Store current warning state
warning('off', 'all');           % Temporarily suppress all warnings
if verLessThan('matlab','9.13')  % MATLAB R2022b is version 9.13
    volshow(RI_3D_temp, 'BackgroundColor', [0 0 0], 'Renderer','MaximumIntensityProjection'); % For MATLAB 2022a and earlier
else
    images.compatibility.volshow.R2022a.volshow(RI_3D_temp, 'BackgroundColor', [0 0 0], 'Renderer','MaximumIntensityProjection');% For MATLAB 2022b and later
end
warning(warnState); 
frame = getframe(fig);
close(fig);

figure('Name','Orthogonal Views','Position', screenSize);
subplot(221); imshow(RI_3D_temp(:,:,zcR),[]); axis equal; axis tight; colorbar; title('x-y view');
subplot(222); imshow(squeeze(RI_3D_temp(:,ycR,:)),[]); axis equal; axis tight; colorbar; title('x-z view');
subplot(223); imshow(squeeze(RI_3D_temp(xcR,:,:)).',[]); axis equal; axis tight; colorbar; title('y-z view');
subplot(224); imshow(frame.cdata,[]); axis equal; axis tight; title([title_temp,' RI Volume (MIP)']);

%% -------------------------------------------------------------------------
% Step 5: Save results
% -------------------------------------------------------------------------
if saveResults
    if useAbeCorrection
        CTF_abe = CTF_abe.*(abs(CTF_abe)>10^-2);
    else
        CTF_abe = CTF;
        zernikeCoeff = [];
    end
    
    saveDir = 'Result_simulation';
    if ~exist(saveDir,'dir')
       mkdir(saveDir);
    end
   
    if useDarkfield
        file_name = fullfile(saveDir,[char(sampleName),'_AFP_simulation_',char(approxMethod),'_','Darkfield','.mat']);
        save(file_name,'RIGT_pad_LP_fullfield','RI_3D_fullfield','RIGT_pad_LP','RI_3D','CTF_abe','CTF_GT','zernikeCoeff','zernikeGT','-v7.3');
    else
        file_name = fullfile(saveDir,[char(sampleName),'_AFP_simulation_',char(approxMethod),'_','NAmatching','.mat']);
        save(file_name,'RIGT_pad_LP','RI_3D','CTF_abe','CTF_GT','zernikeCoeff','zernikeGT','-v7.3');
    end
end

