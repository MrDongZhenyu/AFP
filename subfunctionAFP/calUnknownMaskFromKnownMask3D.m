function [unknownMask,linearArea,knownFT,linearIdx,linearIdx_full] = calUnknownMaskFromKnownMask3D(recMask,fullMask,ftRecons,k_illu,k_carrier,kzNormMap,indexingCTF,X,Y,xc,xcR,yc,ycR,zcR,xres,zres,varargin)
% This function calculates the unknown part position of the entire 
% <spectrum> using the known part, the third input can be the indicator of 
% what <spectrum> is measured.
%
% Notes: 1. The <spectrum> can be the 3D Fourier spectrum of total field or 
%        scattering potential. 
%        2. Since the 3D-CTF has only one z-value for each (x, y) position, 
%        it can be treated as the projection along the z-axis, effectively 
%        making it a 2D spectrum. We simply project it back to 3D when
%        needed.
%   
% Inputs:
%   1. recMask: A 3D binary mask, represents the prior spectrum region 
%   2. fullMask: A 2D binary mask for CTF 
%   3. ftRecons: 3D complex-valued spectrum, represents the prior spectrum
% Inputs (varargin):
%   1. mip: use maximum intensity projection to determine the linear area 
%           in the reconstrution
%   2. darkfield: Whether image belongs to darkfield image
%   3. tol: For the FST prior. When the ratio of (the amplitude of) the 
%           crossterm over the other parts is below this tolerence, 
%           it is approximated as the linear region. 
%   4. thickness: The sample thickness in pixel, used for the FST prior.
%   5. use_gpu: Whether to use gpu for the calculation
%
% Outputs:
%   1. unknownMask: A 2D binary mask, represents the unknown position of 
%                   the <spectrum> in its 2D projection.
%   2. linearArea: A 2D binary mask for the linear region in the Fourier 
%                  spectrum of the intensity measurement.
%   3. knownFT: The known part of the <spectrum> (2D)
%   4. linearIdx: The linear index (position) of the unknown part of the 
%                 <spectrum> in a 3D matrix of size(recMask), when this
%                 <spectrum> is projected back into 3D frequency space. 
%   5. linearIdx_full: The linear index (position) of the entire <spectrum> 
%                 (the Ewald's sphere) in a 3D matrix of size(recMask), when 
%                 this <spectrum> is projected back into 3D frequency space.
%
% By Ruizhi Cao, Zhenyu Dong

% Default parameters
useMIP = false;     % whether to use maximum intensity projection
isDark = false;     % whether image belongs to darkfield image
tol    = 0.1;       % It will be updated if it is specified in the input
usethickness = false;

if ~isempty(varargin)
    idx = 1;
    while idx <= length(varargin)
        switch lower(varargin{idx})
            case {'max projection','max amp','mip'}
                maxAmpMap = varargin{idx+1};
                useMIP = true;
                idx = idx+2;
            case {'tol','tolerence'}
                tol = varargin{idx+1};
                idx = idx+2;
            case {'gpu','use_gpu'}
                use_gpu = varargin{idx+1};
                idx = idx+2;
            case {'sample thickness','thickness'}
                thickness = varargin{idx+1};
                usethickness = true;
                idx = idx+2;
            case {'darkfield','dk','dark field'}
                isDark = true;
                idx = idx+1;
            otherwise
                error(['Supported options: ''max projection'', ''use_gpu''' ...
                    '''tolerence'',''thickness'',''darkfield''.']);
        end
    end
end

k2use = round(k_illu);

%% Incorporate the FST prior
if usethickness && thickness < size(fullMask,1)
    % The window function is rectangular shape, with a width equal to 
    % the sample thickness.
    zSize = length(recMask(1,1,:));
    windowZ = zeros(zSize*size(indexingCTF,1)/size(recMask,1),1);
    windowZ(1:thickness) = 1;

    % Ideally, this should be squared to get the convolution in the Fourier space.
    % However, it is the same for this specific window function
    ftWindow = fft(windowZ);
    ftWindow = abs(ftWindow)/ftWindow(1);
    pixelTolZ = find(ftWindow<(1-tol),1)-1;

    addMask = false(size(recMask));
    addFT = zeros(size(recMask));

    if use_gpu
        recMask = gpuArray(recMask); 
        ftRecons = gpuArray(ftRecons);
        addMask = gpuArray(addMask);
        addFT = gpuArray(addFT);
    end

    diffTemp = diff(recMask>10^-3,1,3); % calculate the 'gradient' in z
    diffNeg = padarray(diffTemp == -1,[0,0,1],false,'post'); % The mask lies mostly at the Ewald's sphere surface
    diffPos = padarray(diffTemp == 1, [0,0,1],false,'pre');  % The mask lies mostly at the central z-plane 
    maskNeg = (recMask>10^-3).*diffNeg;
    maskPos = (recMask>10^-3).*diffPos;
    
    % Mask out the spectrum
    shellNeg = ftRecons.*diffNeg;
    shellPos = ftRecons.*diffPos;
    
    % We aim to expand Ewald's sphere in 3D by assuming that the spectrum 
    % remains the same when it is shifted by up to (pixelTolZ-1) pixels 
    % along the z-frequency axis
    for idx = 1:pixelTolZ-1
        % Positive (xxxPos) part (the third term) considers brightfield spectrum extension. 
        % It has no influence on darkfield spectrum extension
        addMask = addMask | circshift(maskNeg,[0,0,idx]) | circshift(maskPos,[0,0,-idx]); 
        addFT = addFT + circshift(shellNeg,[0,0,idx]) + circshift(shellPos,[0,0,-idx]);
    end
    
    % Expand the spectrum and the mask
    ftRecons = ftRecons + addFT.*(1-recMask>10^-3);
    recMask = (addMask & ~recMask) | recMask;
end

%% Find the linear part within the intensity measurement's Fourier spectrum
kzOffset = (k_carrier - sqrt(k_carrier^2 - (k_illu(1)*xres)^2 - (k_illu(2)*xres)^2))/zres;

kzRaw = kzNormMap + kzOffset + zcR;
zGridLib = (1:zSize).';
kzRecBd = squeeze(max(permute(maskNeg,[3,1,2]).*zGridLib,[],1)) + pixelTolZ-1;
bd = calBoundary([xcR-k2use(1),ycR-k2use(2)],size(kzRaw));
temp = kzRaw - kzRecBd(bd(1):bd(2),bd(3):bd(4)).*fullMask;
knownMask = (temp < 0.6) & fullMask;
unknownMask = fullMask & ~knownMask;

if use_gpu
    unknownMask = gather(unknownMask);
    kzRecBd = gather(kzRecBd);
end
 
% Remove small regions that cause discontinuity
if ~isempty( round(0.1*max([regionprops(unknownMask,'Area').Area]))) 
    unknownMask = bwareaopen(unknownMask, round(0.1*max([regionprops(unknownMask,'Area').Area])));
end
knownMask = fullMask & ~unknownMask & (kzRecBd(bd(1):bd(2),bd(3):bd(4))>0);
unknownMask = fullMask & ~knownMask;

if isDark
    disToLine = ((X-xc).*k2use(1) + (Y-yc).*k2use(2))/norm(k2use);
    slit = abs(disToLine) < 2;
    if sum(sum(unknownMask.*slit)) < 5 % 5 pixel tolerence
        unknownMask = unknownMask & (disToLine<0);
        knownMask = fullMask & ~unknownMask;
    end
end

% position of the entire Ewald's sphere
kz2use = round(kzNormMap + kzOffset) + zcR;
linearIdx_full =  sub2ind(size(recMask),X(indexingCTF)-k2use(1) - xc + xcR,...
                                        Y(indexingCTF)-k2use(2) - yc + ycR,...
                                        kz2use(indexingCTF));

% position of the unknown part of the Ewald's sphere
kz2use = min(kzRecBd(bd(1):bd(2),bd(3):bd(4)).*knownMask,kz2use.*knownMask) + kz2use.*(~knownMask);
linearIdx = sub2ind(size(recMask),X(indexingCTF)-k2use(1) - xc + xcR,...
                                  Y(indexingCTF)-k2use(2) - yc + ycR,...
                                  kz2use(indexingCTF));
knownFT = zeros(size(fullMask));
knownFT(indexingCTF) = ftRecons(linearIdx);
knownFT = knownFT.*knownMask;

if useMIP
    bd2use = calBoundary([xcR - k2use(1),ycR - k2use(2)],size(fullMask));
    tempMap = maxAmpMap(bd2use(1):bd2use(2),bd2use(3):bd2use(4));
    temp = (tempMap ~= 0) & unknownMask;
    tempMap = ((tempMap == 0) & unknownMask)*sum(sum(tempMap.*temp))./sum(temp(:))+ tempMap;
else
    tempMap = true(size(knownMask));
end

% Find the linear region
unknownCorrTemp = fftshift(fft2(unknownMask.*tempMap));
unknownCorr = fftshift(ifft2(ifftshift(unknownCorrTemp.*conj(unknownCorrTemp))));

crossCorrTemp = fftshift(fft2(knownMask.*tempMap));
crossCorr = fftshift(ifft2(ifftshift(unknownCorrTemp.*conj(crossCorrTemp))));

linearArea = crossCorr > 1/tol*(unknownCorr + 0.001 + circshift(rot90(crossCorr,2),mod(size(crossCorr)+1,2)));

end

