function [idxOverlapRef, idxOverlapComp,varargout] = calOverlap3D(recFTframe,mycoord,nPairs2use_eachSpectrum,maxCTF,idxMeas2use,pixelSizeZ,kzMap,varargin)
% This function is used to calculate CTF overlap in 3D.
% Theoretically, two Ewald's sphere (3D-CTF) only have overlap along an arc
% shape in 3D, correspond to a line shape in the 2D projection along z axis 
% (we call it: exact overlap line later). Here, we incorporate the concept
% of Finite Sample Thickness (FST) prior to enlarge this overlapping region, 
% which can help to reduce the number of NA-matching measurements required 
% to obtain a good solution. 

% Note: Variables named kzxxx represent the wavenumber values calculated 
%       in the z-frequency space

% Inputs (varargin):
%   1. edgePixel: The number of pixels we want to exclude from the exact 
%                 overlap line (halfed/one-sided).
%   2. thickness: Sample thickness in pixels, related to FST prior.
%   3. similarity tolerence: If kz difference between two points slightly 
%                 staggered on the Ref and Comp Ewald's sphere is smaller 
%                 than a Tol value, it can also be approximated as overlap.
%   4. lowFreqWeights: We already calculated this in the main aberration 
%                 correction function. But we want to extract the pixel 
%                 values only within the overlap region.
%   5. useAdaptiveRegion: Whether to use FST prior or not. Set to 'True' to
%                         enable FST prior.

% Outputs:
%   1. idxOverlapRef: The relative (x,y) coordinate of the pixels (compared 
%                     to (xc,yc)) in the overlap of the reference image
%   2. idxOverlapComp: ~ of the other image in the pair.
% Optional outputs: (Note that they are ordered as well)
%   1. nPixelPerImg: Number of overlapping pixels for each pair
%   2. nPixelTol: Same as it is defined in the main aberration correction
%                 function. 
%   3. kzDiffLib: The 'kzDiff' value of all the pixels in the overlap region,
%                 see the explanations from line 95-107.
%   4. lowFreqWeights: This is used to generate a weight distribution for
%                      the pixels in the overlap region, based on the 
%                      spatial frequency position of the overlap pixels, 
%                      with lower weight values assigned to pixels with 
%                      lower frequency values.
%  5. locZeroFreq: Location of the zero spatial frequency in the overlap.
%                  It is used primarily for internal test and not useful 
%                  in the final code.

% Requirement: Symbolic Math Toolbox
% By Ruizhi Cao, Zhenyu Dong, Haowen Zhou

useLowFreqWeights = false; % In default, do not apply this weight function
useAdaptiveRegion = true;  % In default, FST prior is used.

if ~isempty(varargin)
    idx = 1;
    while idx <= length(varargin)
        switch lower(varargin{idx})
            case {'pixel tolerence','tolerence'}
                nPixelTol = varargin{idx+1};
                idx = idx+2;
            case {'edge','edge pixel'}
                edgePixel = varargin{idx+1};
                idx = idx+2;
            case {'thickness','sample thickness'}
                usethickness = true;
                thickness = varargin{idx+1};
                idx = idx + 2;
            case {'similarity tolerence'}
                similarityTol = varargin{idx+1};
                idx = idx+2;
            case {'lowfreqweights'}
                lowfreqWeights =  varargin{idx+1}; 
                useLowFreqWeights = true;
                idx = idx+2;
            case {'useadaptiveregion'}
                useAdaptiveRegion =  varargin{idx+1}; 
                idx = idx+2;
            otherwise
                error(['Supported options are: ''pixel tolerence'', ' ...
                    '''thickness'', ''similarity tolerence'' , ' ...
                    '''lowfreqweights'', and ''edge pixel''.']);
        end
    end
end

if ~exist('nPixelTol','var')
    nPixelTol = 1;
end

if ~exist('edgePixel','var') 
    edgePixel = 0;
end

if nPixelTol > 3
    warning('Pixel tolerence is larger than 3, might lead to poor estimation');
end

if usethickness
    % The FST prior is used here to determine the maximum difference in kz 
    % values between two points with the same x-y frequency on two Ewald's 
    % spheres (from two NA-matching measurements). 
    % When the kz difference is smaller than 'kzDiffTol', it can be assumed 
    % that the sample phase terms are identical between the two points. 
    % Their (x,y) coordinates are recorded and their phase value difference 
    % will be used for aberration correction.

    % 'kzDiffTol' is determined by finding the half-width of a sinc function 
    % that exceeds the 'similarityTol' value. This sinc function is the 
    % Fourier transform of a rectangular function, with the rectangular width 
    % corresponding to the sample thickness (detailed in Supplementary Note 1).

    syms kzDiffTolx
    eqn = sinc(thickness*pixelSizeZ*kzDiffTolx) == similarityTol;
    searchRange = [0, Inf];
    kzDiffTol = double(vpasolve(eqn,kzDiffTolx,searchRange));
    clear kzDiffTolx;
end

[~,~,nBrightField] = size(recFTframe);
% indices for spectrums that might have overlap with the selected reconstructed spectrum
idxLib = mod(idxMeas2use+(1:nPairs2use_eachSpectrum)-1,nBrightField)+1; 

idxOverlapTemp1 = cell(1,length(idxLib));
idxOverlapTemp2 = cell(1,length(idxLib));
kzDiffLibTemp = cell(1,length(idxLib));
if useLowFreqWeights
    lowfreqMaskLib = cell(1,length(idxLib));
end
nPixelPerImg = zeros(1,length(idxLib)); % number of overlapping pixels for each pair

idxCount = 1;
locZeroFreq = zeros(length(idxLib),2);  % the zero frequency (DC) position for each sub-spectrum 
for idxNext = idxLib
    % Get the (x,y) coordinates for pixels on the exact overlap line 
    % between two Ewald's spheres, (0,0) in the image center.
    relShift = mycoord(idxMeas2use,:) - mycoord(idxNext,:);
    normRelShift = norm(relShift);
    lengthOverlapLine = max(sqrt(maxCTF.^2 - (normRelShift/2).^2) - edgePixel, 1);  % half of the total length
    directionLine = [relShift(2),-relShift(1)]/normRelShift;
    offsetLine = (mycoord(idxMeas2use,:) + mycoord(idxNext,:))/2;
    lineIdxTemp = -round(lengthOverlapLine):1:round(lengthOverlapLine);
    idx2D = lineIdxTemp.'*directionLine - offsetLine;  
    
    if useAdaptiveRegion
        % FST prior is used to define overlap region
        kzRef = circshift(kzMap,-mycoord(idxMeas2use,:));
        kzComp = circshift(kzMap,-mycoord(idxNext,:));
        kzDiff = abs(kzRef - kzComp);
        overlapRegion = (kzDiff < kzDiffTol).*((kzRef > 0)|(kzComp > 0));
    else
        % Only use pixels on the exact overlap line as overlap region
        overlapRegion = (kzDiff == 0);
    end
    
    if useLowFreqWeights
        lowfreqMaskLib{idxCount} = lowfreqWeights(overlapRegion == 1);
    end
    
    % Find the (x,y) coordinates for each pixel in the overlap region. 
    % These coordinates have different values for the two Ewald's spheres, 
    % resulting in 'idxOverlapTemp1' and 'idxOverlapTemp2'. We also store 
    % the 'kzDiff' value of all the pixels in the overlap region 
    % (named as 'kzDiffLibTemp').
    [idx1DxPerImg,idx1DyPerImg] = find(overlapRegion == 1);
    sizeOvp = size(overlapRegion);
    idx2DPerImg = [idx1DxPerImg - floor(sizeOvp(1)/2+1), idx1DyPerImg  - floor(sizeOvp(2)/2+1)];
    idxOverlapTemp1{idxCount} = [round(idx2DPerImg + mycoord(idxMeas2use,:)),...
                                 idxMeas2use*ones(size(idx2DPerImg,1),1)];
    idxOverlapTemp2{idxCount} = [round(idx2DPerImg + mycoord(idxNext,:)),...
                                 idxNext*ones(size(idx2DPerImg,1),1)];
    kzDiffLibTemp{idxCount} = kzDiff(overlapRegion == 1);
    nPixelPerImg(idxCount) = size(idx2DPerImg,1);

    if nargout == 6
        absVal = vecnorm((idxOverlapTemp1{idxCount}(:,1:2)).' - mycoord(idxMeas2use,:).');
        [minVal,minIdx] = min(absVal);
        if minVal > 2.5^2
            error(['Zero frequency is outside the intersection. Currect idx = ',num2str(idxMeas2use),', idxNext = ',num2str(idxNext)]);
        else
            locZeroFreq(idxCount,:) = [rem(minIdx-1,size(idx2D,1))+1,floor((minIdx-1)/size(idx2D,1))+1];
        end
    end

    idxCount = idxCount + 1;  
end

% For each output variable, concatenate all pairs into a single matrix
nPixelOverlap = sum(nPixelPerImg);
idxOverlapRef   = zeros(nPixelOverlap,3);
idxOverlapComp  = zeros(nPixelOverlap,3);
kzDiffLib = zeros(nPixelOverlap,1);
if useLowFreqWeights
    lowFreqWeights = zeros(nPixelOverlap,1);
end
idxSt = 1;
for idxOverlap = 1:length(idxLib)
    idxOverlapRef(idxSt:idxSt + nPixelPerImg(idxOverlap) - 1,:) = idxOverlapTemp1{idxOverlap};
    idxOverlapComp(idxSt:idxSt + nPixelPerImg(idxOverlap) - 1,:) = idxOverlapTemp2{idxOverlap};
    kzDiffLib(idxSt:idxSt + nPixelPerImg(idxOverlap) - 1,1) = kzDiffLibTemp{idxOverlap};
    if useLowFreqWeights
        lowFreqWeights(idxSt:idxSt + nPixelPerImg(idxOverlap) - 1,1) = lowfreqMaskLib{idxOverlap};
    end
    idxSt = idxSt + nPixelPerImg(idxOverlap);
end

switch nargout
    case 3
        varargout{1} = nPixelPerImg;
    case 4
        varargout{1} = nPixelPerImg;
        varargout{2} = nPixelTol;
    case 5
        varargout{1} = nPixelPerImg;
        varargout{2} = nPixelTol;
        varargout{3} = kzDiffLib;
    case 6
        varargout{1} = nPixelPerImg;
        varargout{2} = nPixelTol;
        varargout{3} = kzDiffLib;
        varargout{4} = lowFreqWeights;
    case 7
        varargout{1} = nPixelPerImg;
        varargout{2} = nPixelTol;
        varargout{3} = kzDiffLib;
        varargout{4} = lowFreqWeights;
        varargout{5} = locZeroFreq;
end

end

