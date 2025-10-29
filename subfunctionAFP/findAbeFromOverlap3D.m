function [CTF_abe,zernikeCoeff] = findAbeFromOverlap3D(recFTframe,myShift,pixelSizeZ,CTF,kzMap,varargin)
% This function estimates the aberration from multiple KK reconstructions 
% with different NA-matching illumination angles
% 
% Inputs:
%   1. recFTframe: The Fourier spectra from different NA-matching 
%                  illuminations
%   2. myShift: The coordinates of illumination angles in pixels  
%   3. kzMap: unit: 1/um. z-axis wavenumber for each point in the Fourier 
%             domain, it has the same shape as the Ewald's sphere. 
% 
% Input Options:
%   1. zernike mode: Zernike mode to optimize for aberration estimation
%   2. use image: Select which image to use in the algorithm. By default, 
%                 all input images are used.
%   3. weighted: Use a weighted matrix in the algorithm to rank the
%                pixels in the overlap region, allowing the algorithm to 
%                prioritize stronger and more reliable signals. 
%                We have four aspects of consideration for the weights.
%                See lines 299-329 below for details.
%   4. max CTF: Define the maximal spatial frequency of the CTF manually. 
%               Use this for a more accurate (subpixel) cutoff frequency 
%               estimation.
%   5. closest n pairs: For a given spectrum, pair it with the closest n 
%                       other spectrums and extract their overlaps for 
%                       aberration correction. By default, the program 
%                       finds the closest N/3 pairs, where N is the number 
%                       of measurements. Change variable 
%                       "nPairs2use_eachSpectrum" for the defualt #.
%   6. regularizer: Regularization term, no need to tune, very minor 
%                   influence to the final result.
%   7. timer: Show computation time remaining.
%   8. thickness: Sample thickness in pixels, will be used for FST prior in
%                 'calOverlap3D.m' function;
%   9. tolerence: For the FST prior calculation, as detailed in the 
%                 'calOverlap3D.m' function, we use a default value of 
%                 sqrt(2)/2. This tolerence value determines the half-width 
%                 of a function (such as gaussian or sinc) when it drops to 
%                 the 'tolerence' percentage of its peak value.
%  10. vis_samplingmap: Whether to visualize the aberration sampling map.
%                 (The sampling density distribution for pixels in the CTF)
%
% Outputs:
%   1. CTF_abe: Complex-valued aberration function
%   2. zernikeCoeff: Zernike coeffients starting from the piston term
%
% By Ruizhi Cao, Zhenyu Dong, Haowen Zhou

%% Set parameters and preparation steps
% Default parameters
[xsize,ysize,nBrightField] = size(recFTframe);
useimg   = 1:nBrightField;
nPairs2use_eachSpectrum = floor(nBrightField/3);  % Define the default number of pairs for each measurment.
zernikeMode2Opt = 3:35;         % Zernike modes to be corrected
regZernike = 0.001;
marginPix = 2;
nPixelTol = 0;                  % Tolerence for the spectrum overlap. A larger tolerence means 
                                % that a larger area is used to determine the aberration.
usetimer = true;
useWeights = true;
desiredMed2MaxWeightRatio = 2;  % Ratio of weight of the median signal to that of the maximal signal 
mycoord = round(myShift);
CTFThreshold = 1*10^-2;
similarityTol = sqrt(2)/2;
vis_samplingmap = false;      

% load input parameters
idx = 1;
while idx<=length(varargin)
    switch lower(varargin{idx})
        case {'zernikemode','zernike mode'}
            zernikeMode2Opt = varargin{idx+1};
            idx = idx + 2;
        case {'useimage','use image'}
            useimg = varargin{idx+1};
            idx = idx + 2;
        case {'use weight','weighted average','weighted','weight'}
            useWeights = varargin{idx+1};
            idx = idx + 2;
        case {'max ctf','cutoff'}
            maxCTFUser = varargin{idx+1};
            idx = idx + 2;
        case {'find overlap','find pairs','closest n pairs','n pairs'}
            if isnumeric(varargin{idx+1}) && varargin{idx+1}>0
                nPairs2use_eachSpectrum = varargin{idx+1};
                idx = idx + 2;
            else
                error(['Please specify the argument for ''closest n ' ...
                       'pairs'' with a natural number.']);
            end
        case {'regularizer','reg'}
            regZernike = varargin{idx+1};
            idx = idx + 2;
        case {'timer','use timer','usetimer'}
            usetimer = varargin{idx+1};
            idx = idx + 2;
        case {'thickness','sample thickness'}
            SamThick = varargin{idx+1};
            idx = idx + 2;
        case {'similarity tolerence','tol','tolerence'}
            similarityTol = varargin{idx+1};
            idx = idx+2;
        case {'vis_samplingmap'}
            vis_samplingmap = varargin{idx+1};
            idx = idx+2;
        otherwise
            error(['Supported options are ''use image'', ''weighted'', ' ...
                '''cutoff'', ''closest n pairs'', ''regularizer'', ' ...
                '''usetimer'', ''thickness'', ''similarity tolerence'', ' ...
                '''vis_samplingmap'', and ''zernike mode''.']);
    end
end

% Check for valid zernike modes
if any(zernikeMode2Opt<3)
    warning('The first 3 zernike mode is dropped, as those produce a worse estimation.');
    zernikeMode2Opt = zernikeMode2Opt(zernikeMode2Opt>=3);
end

% Validate useimg
if any(useimg<1) || any(useimg>nBrightField)
    error('The index of image must be positive, and should not be greater than the number of images.');
end
useimg = sort(useimg);

% Check nPairs2use
if nPairs2use_eachSpectrum > round(nBrightField/2)
    warning('The program is asked to pair one spectrum with more than half of the acquried spectrums. This is likely to be wrong');
end

% Calculate the center coordinates
xc = floor(xsize/2+1); 
yc = floor(ysize/2+1);
nZernike = length(zernikeMode2Opt);
maxCTF = find(abs(CTF(xc,yc+1:end))<CTFThreshold,1);

% Generate a weight function according to the spatial frequency value
[Y,X] = meshgrid(1:ysize,1:xsize);
R = abs(X-xc + 1i*(Y-yc));
rWeights = log(R*18/20/(maxCTF*0.5)+1)/5.889; % first 5 percent goes 1/2; 5.889 = log(18*20+1);

% Calculate dcAmp
dcAmp = zeros(1,nBrightField);
for idx = 1:nBrightField
    dcAmp(idx) = abs(recFTframe(xc+mycoord(idx,1),yc+mycoord(idx,2),idx));
end

% Find the maximum DC amplitude
maxDC = max(dcAmp);

% Determine the resolution pixel size
resPix = min([xsize,ysize]) - 1 - 2*maxCTF;

% Adjust margin pixel size if necessary
if marginPix*2>resPix
    marginPix = floor(resPix/2);
end

% Check for existence of maxCTFUser and compare with maxCTF
if exist('maxCTFUser','var')
    if maxCTF < maxCTFUser
        warning(['The thresholding-based cutoff frequency estimate is ' ...
                 'smaller than the given value. The code will overwrite ' ...
                 'this value. Please check the given cutoff freq.']);
        maxCTFUser = maxCTF;
    end
    temp = linspace(-maxCTF/maxCTFUser,maxCTF/maxCTFUser,2*maxCTF+1+2*nPixelTol);
else
    temp = linspace(-1,1,2*maxCTF+1+2*nPixelTol);
end

% Calculate the Zernike polynomial matrix (Hz)
[Yz,Xz] = meshgrid(temp,temp);
[theta,r] = cart2pol(Yz,Xz);
idx2use = (r<=1);
zernikeTemp = zernfun2(zernikeMode2Opt,r(idx2use),theta(idx2use));
Hz = zeros((2*maxCTF+1+2*nPixelTol)^2,nZernike); % zernike operator that generates aberration
for idx = 1:nZernike
    Hz(idx2use,:) = zernikeTemp;
end
clear zernikeTemp theta r;

%% Build matrix A for the linear equation y = Ax
% Allocate placeholder for variables
Hdiff = [];        % operator that calculates the aberration differences
                   % note here Hdiff is basically D_{il} in our derivation.
numPhaseMeas = []; % number of phase measuerement in each pair
phaseMeas = [];    % phase difference values
ncolHdiff = (2*maxCTF+1+2*nPixelTol)^2; % number of columns for matrix 'Hdiff'
weightsVct = [];   % the total weights matrix
offsetIdx = [];    % index of phase offset term location
weightsKzDiff = []; % weights related to 'kzDiff'
weightsLowFreq = []; % weights related to position in spatial frequency

if vis_samplingmap
    overlapIdxRef_total = [];
    overlapIdxComp_total = [];
end

if usetimer
    myTimer = CalTimeRemain([length(useimg),nPairs2use_eachSpectrum]);
end

% Build linear equation
for idx = useimg
    % (xc,yc) is the (0,0) position in overlapIdxRef and overlapIdxComp
    [overlapIdxRef,overlapIdxComp,nPixelPerImg,nPixelTol,kzDiffLib,lowFreqTemp] = calOverlap3D(recFTframe,mycoord,nPairs2use_eachSpectrum,maxCTF,idx,pixelSizeZ,kzMap,...
                                                                                               'tolerence',nPixelTol,'thickness',SamThick,'useadaptiveregion',true,...
                                                                                               'lowfreqweights',rWeights,'similarity tolerence',similarityTol); 
    linearIdxRef = sub2ind([xsize,ysize,nBrightField], overlapIdxRef(:,1)  + xc,   overlapIdxRef(:,2) + yc,  overlapIdxRef(:,3));
    linearIdxComp = sub2ind([xsize,ysize,nBrightField],overlapIdxComp(:,1) + xc, overlapIdxComp(:,2)  + yc,  overlapIdxComp(:,3));
    
    if useWeights
        weightsKzDiff = [weightsKzDiff;kzDiffLib];
        weightsLowFreq = [weightsLowFreq;lowFreqTemp];
    end
    
    idxEd = cumsum(nPixelPerImg);
    idxSt = [0,idxEd(1:end-1)]+1;
    for idxPair = 1:length(nPixelPerImg)
        idxNext = overlapIdxComp(idxSt(idxPair),3);
        
        % index of the zero-freq in the vectorized spectrum 
        % (on cropped pupil coordinates for Zernike)
        originSpectrum1 = (maxCTF+mycoord(idx,2))*(2*maxCTF+1+2*nPixelTol)     + maxCTF + nPixelTol + mycoord(idx,1) + 1;
        originSpectrum2 = (maxCTF+mycoord(idxNext,2))*(2*maxCTF+1+2*nPixelTol) + maxCTF + nPixelTol + mycoord(idxNext,1) + 1;
   
        posIdx = sub2ind((2*maxCTF+1+2*nPixelTol)*[1,1], overlapIdxRef(idxSt(idxPair):idxEd(idxPair),1)  + maxCTF+1+nPixelTol,...
                                                         overlapIdxRef(idxSt(idxPair):idxEd(idxPair),2)  + maxCTF+1+nPixelTol);
        negIdx = sub2ind((2*maxCTF+1+2*nPixelTol)*[1,1], overlapIdxComp(idxSt(idxPair):idxEd(idxPair),1) + maxCTF+1+nPixelTol,...
                                                         overlapIdxComp(idxSt(idxPair):idxEd(idxPair),2) + maxCTF+1+nPixelTol);
        nMeasure = nPixelPerImg(idxPair);
        
        % the phase of phaseTemp is the total phase difference between 
        % two paired regions 
        phaseTemp = recFTframe(linearIdxRef(idxSt(idxPair):idxEd(idxPair))).*...
                    conj(recFTframe(linearIdxComp(idxSt(idxPair):idxEd(idxPair))));
        
        % get zero spatial frequency position in the linearIdxRef list
        zeroFreq = sub2ind([xsize,ysize],xc+mycoord(idx,1),yc+mycoord(idx,2)) + xsize*ysize*(idx-1);
        zeroFreqInLinearIdx = find(linearIdxRef(idxSt(idxPair):idxEd(idxPair)) == zeroFreq);
        if ~isempty(zeroFreqInLinearIdx)
            zeroFreqList(length(numPhaseMeas)+1) = zeroFreqInLinearIdx + sum(numPhaseMeas,'all');
        end
        
        % We assign weights based on the SNR of the signal; the higher the 
        % SNR for the pair of measurements, the greater the weight applied. 
        % This reflects a higher level of trust for that pair in forming 
        % the linear equation.
        if useWeights
            weightsTemp = abs(phaseTemp);
            weightsVct = [weightsVct;weightsTemp.*(weightsTemp>0)];
        end

        % convert to 2D rectangular shape for phase unwrapping
        phaseTempRect = zeros(xsize,ysize);
        [PutBackIdxX, PutBackIdxY] = ind2sub([xsize,ysize],linearIdxRef(idxSt(idxPair):idxEd(idxPair)) - xsize*ysize*(idx-1));
        PutBackIdxX = PutBackIdxX + mycoord(idx,1);
        PutBackIdxY = PutBackIdxY + mycoord(idx,2);
        PutBackIdx = sub2ind(size(phaseTempRect),PutBackIdxX,PutBackIdxY);
        phaseTempRect(PutBackIdx) = phaseTemp + eps;
        [idx1Dx,idx1Dy] = find( phaseTempRect ~= 0);
        phaseTemp = phaseTempRect(min(idx1Dx):max(idx1Dx),min(idx1Dy):max(idx1Dy));  
        
        % unwrap the phase before aberration extraction
        phaseUnwrapTemp = unwrap2(angle(phaseTemp));
        % phaseUnwrapTemp = unwrap(angle(phaseTemp));  % If unwrap2 function not working, use this line instead. 
        phaseRaw = angle(phaseTemp);

        % force phase unwrap using multiplies of 2pi
        wrappedPhase = phaseUnwrapTemp - phaseRaw;
        x2use = -(2*pi):pi/8:(2*pi);
        N = histcounts(wrappedPhase(:),x2use);
        idx2use = (x2use >= -(pi+pi/4)) & (x2use <= pi+pi/4);
        N(~idx2use) = 0;
        idxPk = find(N == max(N));idxPk = idxPk(1); % use the first peak when multiple maxima are found.
        offsetPk = mean2(wrappedPhase(wrappedPhase>=(x2use(idxPk)-pi/4) & wrappedPhase<=(x2use(idxPk)+pi/4)));
        phaseTemp = phaseRaw + round((wrappedPhase - offsetPk)/(2*pi))*2*pi + offsetPk;

        % extract points within overlapping region (an irregular shape) 
        % from the rectangular shape
        [idx1Dx,idx1Dy] = find(phaseTempRect(min(idx1Dx):max(idx1Dx),min(idx1Dy):max(idx1Dy)) ~= 0);
        MeasureIdx = sub2ind(size(phaseTempRect(min(idx1Dx):max(idx1Dx),min(idx1Dy):max(idx1Dy))), idx1Dx,idx1Dy);

        % prepare matrices for solving the least squares fit
        phaseMeas = [phaseMeas;reshape(phaseTemp(MeasureIdx),[],1)];
        numPhaseMeas = [numPhaseMeas,nMeasure];
        Hdiff = [Hdiff;sparse(repmat(1:nMeasure,[1,2]),[posIdx.',negIdx.'],...
                [ones(1,nMeasure),-ones(1,nMeasure)],nMeasure,ncolHdiff)];
        offsetIdx = [offsetIdx;originSpectrum1,originSpectrum2];

        if usetimer
            myTimer.timeRemain([idx,idxPair]);
        end
    end
    
    if vis_samplingmap 
        overlapIdxRef_total = [overlapIdxRef_total;overlapIdxRef];
        overlapIdxComp_total = [overlapIdxComp_total;overlapIdxComp];
    end
end

% Visualize the aberration sampling map (for debug only)
if vis_samplingmap
    sampling_map = zeros(size(CTF));
    for idx = 1:length(overlapIdxRef_total)
        sampling_map(xc+overlapIdxRef_total(idx,1),yc+overlapIdxRef_total(idx,2)) = sampling_map(xc+overlapIdxRef_total(idx,1),yc+overlapIdxRef_total(idx,2))+1;
        sampling_map(xc+overlapIdxComp_total(idx,1),yc+overlapIdxComp_total(idx,2)) = sampling_map(xc+overlapIdxComp_total(idx,1),yc+overlapIdxComp_total(idx,2))+1;
    end
    figure('Name','Aberration sampling map');imshow(sampling_map,[]);
    viscircles([xc,yc],maxCTF+nPixelTol); 
end

if usetimer
    myTimer.delete;
end

% Populate Hoffset
Hoffset = zeros(sum(numPhaseMeas),nZernike); % operator that calculates the phase offset
% Note: Hoffset is basically D^0_{il} in our derivation
for idx = 1:length(numPhaseMeas)
    idxSt = 1+sum(numPhaseMeas(1:idx-1));
    Hoffset(idxSt:idxSt+numPhaseMeas(idx)-1,:) = repmat(Hz(offsetIdx(idx,1),:) - Hz(offsetIdx(idx,2),:),[numPhaseMeas(idx),1]);
end

%% Apply Weights if useWeights is True (This step is crucial!)
if useWeights
    % 1) Set weight at zero frequencies to zero
    if ~isempty(zeroFreqList)
        weightsVct(nonzeros(zeroFreqList)) = 0;
    end

    % 2) Consider SNR of the signal (higher SNR, higher weight)
    medWeights = median(weightsVct);
    ratioMax2Med = (maxDC.^2)/medWeights;
    factor = ceil(desiredMed2MaxWeightRatio*log10(ratioMax2Med)/(desiredMed2MaxWeightRatio-1));
    weightsVct = log10(weightsVct./max(weightsVct(:))*10^factor +1);
    
    % 3) Apply a linear weight function according to 'kzDiff', remember that
    % we used FST prior to expand the overlap region, the further away from
    % the exact overlap line, the smaller weight we use.
    weightsVct = weightsVct.*linear_weight(weightsKzDiff,min(weightsKzDiff),max(weightsKzDiff),0.2);
    
    % 4) Apply weight according to spatial frequency position of the overlap 
    % points. Pixels with higher spatial frequencies are given larger 
    % weights, as the zero-frequency component is prone to phase 
    % ambiguities due to inaccuracies in illumination angle calibration 
    % or rounding issues.
    weightsVct = weightsVct.*weightsLowFreq;
    
    % Here, Hoverall is our final operator encoded with weights W*D
    Hoverall = weightsVct.*(Hdiff*Hz - Hoffset);  
else
    Hoverall = Hdiff*Hz - Hoffset;
end
clear Hdiff Hoffset;

%% Compensate for floating phase offset: (W*) (D*Z*x + H_e*e) = (W*)y; H_error = H_e'*H_e
lengthPhaseMeas = length(numPhaseMeas);
sigOffset = zeros(lengthPhaseMeas,1);
blockLL = zeros(lengthPhaseMeas,nZernike);
weightsVctSq = weightsVct.^2;   % squared weights
weightsSqSum = numPhaseMeas.';
for idx = 1:lengthPhaseMeas
    idxSt = 1+sum(numPhaseMeas(1:idx-1));
    
    if useWeights
        blockLL(idx,:) = sum(weightsVct(idxSt:idxSt+numPhaseMeas(idx)-1).*Hoverall(idxSt:idxSt+numPhaseMeas(idx)-1,:),1);
        weightsSqSum(idx) = sum(weightsVctSq(idxSt:idxSt+numPhaseMeas(idx)-1));
        sigOffset(idx) = sum(weightsVctSq(idxSt:idxSt+numPhaseMeas(idx)-1).*phaseMeas(idxSt:idxSt+numPhaseMeas(idx)-1));
    else
        blockLL(idx,:) = sum(Hoverall(idxSt:idxSt+numPhaseMeas(idx)-1,:),1);
        sigOffset(idx) = sum(phaseMeas(idxSt:idxSt+numPhaseMeas(idx)-1));
    end
end

%% Least-squares fit of the aberration (see Supplementary Note 1 for details)
Herror = diag(weightsSqSum);
Hnew = [Hoverall'*Hoverall,blockLL';blockLL,Herror];
mtxReg = regZernike*mean(abs(Hnew(:)))*diag([ones(1,nZernike),zeros(1,lengthPhaseMeas)]);
if useWeights
    zernikeCoeff_temp = (Hnew + mtxReg)\[(Hoverall'*(weightsVct.*phaseMeas));sigOffset];
else
    zernikeCoeff_temp = (Hnew + mtxReg)\[(Hoverall'*phaseMeas);sigOffset];
end
zernikeCoeff_new = zernikeCoeff_temp(1:nZernike);

%% Generate the extracted aberration
temp = reshape(Hz*zernikeCoeff_new,[2*maxCTF+1+2*nPixelTol, 2*maxCTF+1+2*nPixelTol]);
CTF_abe = double(CTF);
bdC = calBoundary([xc,yc],[2*maxCTF+1+2*nPixelTol,2*maxCTF+1+2*nPixelTol]);
CTF_abe(bdC(1):bdC(2),bdC(3):bdC(4)) = CTF_abe(bdC(1):bdC(2),bdC(3):bdC(4)).*exp(1i*temp);

%% Extract the coefficients
zernikeCoeff = zeros(1,max(zernikeMode2Opt)+1);
zernikeCoeff(zernikeMode2Opt+1) = zernikeCoeff_new;

if marginPix>0
    expandFactor = (marginPix + maxCTF)/maxCTF;
    fitAbe = imresize(phase_unwrapCG(angle(CTF_abe)),[ceil(xsize*expandFactor),ceil(ysize*expandFactor)],'bilinear');
    bd = calBoundary(floor(size(fitAbe)/2+1),[xsize,ysize]);
    fitAbe = fitAbe(bd(1):bd(2),bd(3):bd(4));
    temp = CTF_abe.*(abs(CTF_abe)>0.01) + (abs(CTF_abe).*(abs(CTF_abe)<=0.01) + 10^-5).*exp(1i*fitAbe);
    CTF_abe =  CTF_abe.*(abs(CTF_abe)>0.01) + (abs(CTF_abe).*(abs(CTF_abe)<=0.01) + 10^-5).*exp(1i*medfilt2(phase_unwrapCG(angle(temp))));
end

end

