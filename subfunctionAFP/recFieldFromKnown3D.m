function [ftRecons,maskRecons,recScatteringPotentialFT_fullfield,mask_fullfield] = recFieldFromKnown3D(imStack,kIllu,ftRecons,maskRecons,recScatteringPotentialFT,CTF_abe,k_carrier,kzNormMap,indexingCTF,zcR,xres,zres,kz,pixelSizeZ,varargin)
% This function reconstructs field and scattering potential based on known/ 
% prior Fourier spectrum of the measurement.
%
% Note: Detailed descriptions of the Finite Sample Thickness (FST) prior
%       can be found in 'calUnknownMaskFromKnownMask3D.m' function.
% 
% Inputs: 1. imStack: Meaurements from the experiment
%         2. kIllu: Illumination vector, which corresponds to the shift 
%                   (in pixel) in spatial frequency domain.
%         3. ftRecons: NA-matching part of the total field Fourier spectrum
%         4. maskRecons: NA-matching Fourier spectrum mask with 1 denotes 
%                        the known FT, 0 the unknown
%         5. recScatteringPotentialFT: NA-matching Fourier spectrum of the 
%                        scattering potential
%         6. CTF_abe: The retrieved complex-valued aberration function
% Input options:
%         1. drift: when enabled (set to true), then the program take 
%                   into consideration that the illumination vector might
%                   be inaccurate
%         2. regularization: Specify the weight of the L2 regularizer
%         3. unknown ratio: Only measurements with unknown spectrum that
%                           are larger than this ratio*CTF will be used.
%         4. high freq threshold (threshold): If the unknown spectrum is 
%                       smaller than this ratio, the newly measured 
%                       spectrum will be averaged, and will be added to the
%                       final spectrum in a later time. The final spectrum 
%                       is expanded immediately when the unkonwn spectrum 
%                       is larger than this ratio.
%         5. timer on: Calculate the remaining time to finish.
%         6. sample thickness: Tell the program the thickness of the
%                       sample. As the darkfield spectrums have disjoint
%                       support with reconstructed spectrum, this code use
%                       the sample thickness as the prior information in
%                       expanding the 3D spectrum using darkfield
%                       measurements.
%         7. darkfield: Specify the measurement type. This tells the
%                       program the measurements are darkfield measurements
%
% Outputs: 1. ftRecons: Total field Fourier spectrum (NA-matching + unknown
%                       part of darkfield)
%          2. maskRecons: Binary mask for the Fourier spectrum coverage 
%                       (NA-matching + unknown part of darkfield)
%          3. recScatteringPotentialFT_fullfield: Fourier spectrum of the 
%                       scattering potential (NA-matching + entire darkfield)
%          4. mask_fullfield: Binary mask for the entire Fourier spectrum 
%                       (NA-matching + entire darkfield)
%
% By Ruizhi Cao, Zhenyu Dong
    
    % Default parameters
    enableDrift      = false; % Whether to take angle calibration error in to consideration
    useTimer         = false;    
    pixelTol         = 2;     % Expand the calculated unknown mask by (roughly) 'pixelTol' pixels when drift is enabled
    highFreqTHLD     = 0.3;   % Threshold for the area of unknown spectrum
    unknownRatio     = 0;     % If the ratio of unknown part to the CTF in area is below this threshold, 
                              % the corresponding measurement is not used in the reconstruction.
    valueOnHold     = false;  % Works with nonzero highFreqTHLD
    userDefinedReg  = false;
    userBrightness  = false;
    autoAmpMatch    = false;  % Whether to correct illumination intensity differences
    isDarkfield     = true;
    marginPixel     = 5;
    replaceMag      = true;   % Whether to match the energy (intensity) of the complex field reconstruction with the actual measurement.
    similarityTol   = 0.2;    % As the darkfield measurements have disjoint support in 3D, this number defines the tolerence
                              % when the spectrum at a slightly different k'_z = (kz+\delta kz) is treated to be approximately the same as the spectrum at kz
                              % If given as input, this is used to generate a sample-thickness-dependent variable (for the FST prior) computed within this function.

    recScatteringPotentialFT_fullfield = recScatteringPotentialFT; % initialize fullfield ScatteringPotential using its NA-matching part
    weightMtx_fullfield = double(maskRecons); % initialize fullfield weightmatrix using its NA-matching part
    
    % load input options
    idx = 1;
    while idx <= length(varargin)
        switch lower(varargin{idx})
            case {'drift','driftcorrection','drift correction'}
                enableDrift = varargin{idx+1};
                if idx+2 <= length(varargin) && isnumeric(varargin{idx+2})
                    pixelTol = varargin{idx+2};
                    idx = idx+3;
                else
                    idx = idx+2;
                end
            case {'reg','regularization','regularizer'}
                myreg = varargin{idx+1};
                userDefinedReg = true;
                idx = idx+2;
            case {'approx','approx_method'}
                approx_method = varargin{idx+1};
                idx = idx+2;
            case {'gpu','use_gpu'}
                use_gpu = varargin{idx+1};
                idx = idx+2;
            case {'unknown ratio','unknownratio','minratio','min ratio'}
                unknownRatio = varargin{idx+1};
                idx = idx+2;
            case {'high freq threshold','threshold','thres'}
                highFreqTHLD = varargin{idx+1};
                if highFreqTHLD > 0.5
                    warning('The threshold should not exceed 0.5. Set to 0.5 instread.');
                    highFreqTHLD = 0.5;
                end
                idx = idx+2;
            case {'timeron','timer on','time'}
                useTimer = true;
                idx = idx+1;
            case {'brightness','intensity'}
                myBrightness = varargin{idx+1};
                userBrightness = true;
                idx = idx + 2;
            case {'intensitycorrection','intensity correction','intensity match'}
                autoAmpMatch = varargin{idx+1};
                idx = idx + 2;
            case {'usedataintensity','dataintensity','data intensity','use data intensity'}
                replaceMag = varargin{idx+1};
                idx = idx + 2;
            case {'sample thickness','thickness'}
                thickness = varargin{idx+1};
                idx = idx+2;
            case {'similarity tolerence'}
                similarityTol = varargin{idx+1};
                idx = idx+2;
            case {'darkfield','dk','dark field'}
                isDarkfield = true;
                idx = idx+1;
            otherwise
                error(['Supported arguments are ''drift'', ''reg'', ' ...
                    '''unknown ratio'', ''threshold'', ''brightness'', ' ...
                    '''intensity correction'', ''use data intensity'', ' ...
                    '''thickness'', and ''timer on''.']);
        end
    end
    
    if ~userDefinedReg
        % when there is no regularization factor specified, the L2 regularizer's
        % weight is chosen based on whether drift correction is enabled
        if enableDrift
            myreg  = 1; % factor for L2 norm regularizer
        else
            myreg  = 0.01;
        end
    end
    
    if isDarkfield
        optsIn = {'darkfield'};
    else
        optsIn = {};
    end
    
    % Get CTF
    maxAmpCTF = max(abs(CTF_abe(:)));
    if  maxAmpCTF < 0.99
        warning('CTF is unnormalized, it will be normalized by default of the program.');
        CTF_abe = CTF_abe./maxAmpCTF;
    end
    CTF = abs(CTF_abe) > 5*10^-3;    
    
    % Checking dimensions of the input arrays
    [xsize,ysize,numIm] = size(imStack);
    if numIm ~= length(kIllu(:,1))
        error('Number of image and number of illumination angles mismatch.');
    end
    
    if userBrightness && length(myBrightness) ~= numIm
        error('Number of intensity calibration points disagrees with the number of images.');
    end
    
    if ~exist('thickness','var')
        thickness = xsize;
    end
    
    % Creating meshgrid for X and Y coordinates
    [Y,X] = meshgrid(1:ysize,1:xsize);
    xc = floor(xsize/2+1);
    yc = floor(ysize/2+1);
    R = abs(X-xc + 1i*(Y-yc));
    
    % Handling marginPixel and boundary calculation
    if marginPixel > 0.1*min(xsize,ysize)
        marginPixel = ceil(max(0.08*min(xsize,ysize),2));
    end
    bdCrop = calBoundary([xc,yc],[xsize-2*marginPixel,ysize-2*marginPixel]);
    
    % Center and size of the reconstruction
    [tempX,tempY,~] = size(ftRecons);
    if tempX ~= tempY
        error('The known reconstruction must be a square image.');
    end
    imsizeRecons = tempX;
    xcR = floor(imsizeRecons/2+1);
    ycR = floor(imsizeRecons/2+1);
    
    if ~exist('maxCTF','var')
        maxCTF = round((sum(CTF(xc,:))-1)/2);
        if sum(CTF(xc,:)) ~= sum(CTF(:,yc))
            error('The input image is not a square image. Please use zero padding.');
        end
    end
    
    if enableDrift
       CTF_larger = (R < maxCTF+1);
    end
    
    areaCTF = sum(CTF(:));
    
    % Initialization
    if highFreqTHLD ~= 0
        unknownTHLD = highFreqTHLD*areaCTF; % threshold in terms of the area of the unknown mask
        ftExpanded = zeros(size(ftRecons)); % store temporary reconstructed scattered field spectrum (with repeats)
        maskExpanded = zeros(size(maskRecons)); % number of repeats for the expanded spectrum
    end

    bd = calBoundary([xcR,ycR],[xsize,ysize]);
    
    if useTimer
        myReconstructionTimer = CalTimeRemain(numIm);
    end
    
    for idx = 1:numIm
        flagSubPixel = false;
        bd2use = bd - repmat(kIllu(idx,:),[2,1]);
        if any(mod(bd2use(:),1)~= 0)
            temp = bd2use - round(bd2use);
            subpixelShift = -[temp(1),temp(3)];
            bd2use = round(bd2use);
            flagSubPixel = true;
        end

        % introduce aberration to the FT to match up with real measurement
        if flagSubPixel
            ampCTF = exact_shift(abs(CTF_abe),-subpixelShift,'real');
            aglCTF = exact_shift(phase_unwrapCG(angle(CTF_abe)),-subpixelShift,'real');
            ampCTF = ampCTF.*(ampCTF>1e-3);
            ampCTF = ampCTF.*(ampCTF<=1) + (ampCTF>1);
            CTF2use = ampCTF.*exp(1i*aglCTF);
        else
            if idx == 1
                CTF2use = CTF_abe.*(abs(CTF_abe)>10^-3);
                ampCTF = abs(CTF2use);
            end
        end

        [unknownMask,linearArea,lowFT,linearIndex3D,linearIndex3D_full] = calUnknownMaskFromKnownMask3D(maskRecons,CTF,ftRecons,kIllu(idx,:),k_carrier,...
                                                                                                        kzNormMap,indexingCTF,X,Y,xc,xcR,yc,ycR,zcR, ...
                                                                                                        xres,zres,'MIP',max(abs(ftRecons),[],3), ...
                                                                                                        'thickness',thickness,'gpu',use_gpu,...
                                                                                                        'tol',similarityTol,optsIn{:});
        lowFT = lowFT.*CTF2use;
        
        if highFreqTHLD ~= 0  &&  sum(unknownMask(:)) > unknownTHLD
            % expand the spectrum of the reconstruction
            ftRecons = ftRecons + ftExpanded./(maskExpanded + eps);
            maskRecons = maskRecons + (maskExpanded ~= 0);
            
            % reset temporary spectrum and its weight mask
            maskExpanded(:) = 0;
            ftExpanded(:) = 0;
            valueOnHold = false;
            
            % regenerate the known spectrum, the known and unknown masks
            [unknownMask,linearArea,lowFT,linearIndex3D,linearIndex3D_full] = calUnknownMaskFromKnownMask3D(maskRecons,CTF,ftRecons,kIllu(idx,:),k_carrier,...
                                                                                                            kzNormMap,indexingCTF,X,Y,xc,xcR,yc,ycR,zcR, ...
                                                                                                            xres,zres,'MIP',max(abs(ftRecons),[],3), ...
                                                                                                            'thickness',thickness,'gpu',use_gpu,...
                                                                                                            'tol',similarityTol,optsIn{:});
            lowFT = lowFT.*CTF2use;
        end
        
        % figure;imshow(unknownMask,[]); % debug purpose, visualize masks for unknown spectrum
        if sum(unknownMask(:)) > unknownRatio*areaCTF
            fieldKnown = ifft2(ifftshift(lowFT));
            imKnown = fieldKnown.*conj(fieldKnown);

            % Calculate the correlation to correct the intensity
            if autoAmpMatch
                corrF = sum(sum((imKnown(bdCrop(1):bdCrop(2),bdCrop(3):bdCrop(4)) - mean2(imKnown(bdCrop(1):bdCrop(2),bdCrop(3):bdCrop(4))))...
                                .*(imStack(bdCrop(1):bdCrop(2),bdCrop(3):bdCrop(4),idx) - mean2(imStack(bdCrop(1):bdCrop(2),bdCrop(3):bdCrop(4),idx)))))./ ...
                        sum(sum((imStack(bdCrop(1):bdCrop(2),bdCrop(3):bdCrop(4),idx) - mean2(imStack(bdCrop(1):bdCrop(2),bdCrop(3):bdCrop(4),idx))).^2));
            else
                corrF = 1;
            end
            
            if userBrightness
                corrF = corrF*myBrightness(idx);
            end
            
            % Subtract the known spectrum
            imReal = imStack(:,:,idx)*corrF - imKnown;
            ftImSub = fftshift(fft2(imReal))*xsize*ysize;
            confinedCTFBoxSize = maxCTF; 

            % Vertices of a smaller square (smaller than the full image) that contains the CTF
            boxVert = [xc - confinedCTFBoxSize,xc + confinedCTFBoxSize,...
                       yc - confinedCTFBoxSize,yc + confinedCTFBoxSize]; 
            
            % Check for aliasing
            if any(boxVert<1) || boxVert(2)>xsize || boxVert(4)>ysize
                error('There is aliasing in the captured image, please reduce the CTF size.');
            end

            if enableDrift
                kernelTol = ones(2*pixelTol+1,2*pixelTol+1);
                unknownMaskOriginal = unknownMask;
                temp = imfilter(unknownMask,kernelTol);
                temp = (temp>0.05*max(temp(:))).*CTF_larger;
                unknownMask(boxVert(1):boxVert(2),boxVert(3):boxVert(4)) = ...
                                temp(boxVert(1):boxVert(2),boxVert(3):boxVert(4));
            end

            % Build convolution matrix 
            Hreduced = calConvMtx(rot90(conj(lowFT(boxVert(1):boxVert(2),boxVert(3):boxVert(4))),2),...
                                  linearArea(xc-confinedCTFBoxSize*2:xc+confinedCTFBoxSize*2,...
                                             yc-confinedCTFBoxSize*2:yc+confinedCTFBoxSize*2),...
                                  unknownMask(boxVert(1):boxVert(2),boxVert(3):boxVert(4)));
            Htemp = Hreduced'*Hreduced;
            absMean = mean2(abs(Htemp));

            % Solve least-squares problem
            if enableDrift
                vecWeight = diag(unknownMaskOriginal(unknownMask)+0.00001)/1.00001; 
                recFTvct = (Htemp + absMean*myreg*vecWeight)\(Hreduced'*ftImSub(linearArea(:) == 1));
            else
                recFTvct = (Htemp + absMean*myreg*eye(size(Htemp,1)))\(Hreduced'*ftImSub(linearArea(:) == 1));
            end

            ftTrue = zeros(xsize,ysize);
            ftTrue(unknownMask) = recFTvct;
            if enableDrift
                unknownMask = unknownMaskOriginal;
                ftTrue = ftTrue.*unknownMaskOriginal;
            end
            
            % Replace magnitude if necessary
            if replaceMag % whether to maintain the (pixel-wise) energy
                fieldTemp = exp(1i*angle(ifft2(ifftshift(ftTrue + lowFT)))).*sqrt(imStack(:,:,idx)*corrF);
                ftTemp = fftshift(fft2(fieldTemp));
                ftTrue(unknownMask) = ftTemp(unknownMask);
            else
                ftTemp = ftTrue + lowFT;
            end
            
            % Correct aberration
            % ftTemp is the entire total field Fourier spectrum retrieved by AFP, 
            % ftTrue is the unknown part of the total field Fourier spectrum.
            ftTemp = ftTemp.*conj(CTF2use)./(abs(CTF2use) + 0.005).*(ampCTF>0.05)*1.005; 
            ftTrue = ftTrue.*conj(CTF2use)./(abs(CTF2use) + 0.005).*(ampCTF>0.05)*1.005; 
            
            % *** Note: one way to check whether the program is correct is
            % to visualize the 'ftTemp' varibale for the bead sample (12 
            % micron diameter for instance) and see whether we can observe 
            % the ring-shaped structures. 
            % figure;imshow(log10(1+abs(ftTemp)),[]); % For debug only.

            % Expand the scattered field spectrum and its mask
            if highFreqTHLD == 0
                % Update ftRecons and maskRecons if highFreqTHLD is zero
                ftRecons(linearIndex3D) = ftRecons(linearIndex3D) + ftTrue(indexingCTF);

                maskRecons(linearIndex3D) = maskRecons(linearIndex3D) + unknownMask(indexingCTF).*(ampCTF(indexingCTF)>0.05);
            else
                % Store the newly reconstucted scattered field spectrum to ftExpanded
                ftExpanded(linearIndex3D) = ftExpanded(linearIndex3D) + ftTrue(indexingCTF);
                
                maskExpanded(linearIndex3D) = ...
                    maskExpanded(linearIndex3D) + unknownMask(indexingCTF).*(ampCTF(indexingCTF)>0.05);
                valueOnHold = true;
            end
            
            % Stitch scattering potential
            if strcmp(approx_method, 'Born')
                % The first order scattered field equals to the scattered field and thus equals to the total field 
                ftTemp = ftTemp;

            elseif strcmp(approx_method, 'Rytov')
                % Also need to pad ftTemp before adding DC term if darkfield extension goes beyond 2*NA_obj
                ftTemp(xc + round(kIllu(idx,1)),yc + round(kIllu(idx,2))) = ftTemp(xc + round(kIllu(idx,1)),yc + round(kIllu(idx,2))) + ftRecons(xcR,ycR,zcR);
                ftTemp = padarray(ftTemp,[xsize/2,ysize/2],0,'both');
                pad_factor = numel(ftTemp)/xsize/ysize;
                ftTemp = ftTemp*pad_factor;
                temp = ifft2(ifftshift(circshift(ftTemp,-round(kIllu(idx,:)))));
                phs = unwrap2(angle(temp));
                amp = abs(temp);
                ftTemp = circshift(fftshift(fft2((log(amp)+1j*phs))),round(kIllu(idx,:)))/pad_factor;
                ftTemp = ftTemp(xsize/2+1:end-xsize/2,ysize/2+1:end-ysize/2).*((ampCTF.*indexingCTF)>0.05);
            end
            
            % Convert to scattering potential using Fourier diffraction theorem
            ftTemp = ftTemp*-2j.*kz*(2*pi*(size(maskRecons,1)/size(CTF,1))^3/pixelSizeZ);
            recScatteringPotentialFT_fullfield(linearIndex3D_full) = recScatteringPotentialFT_fullfield(linearIndex3D_full) + ftTemp(indexingCTF);
            weightMtx_fullfield(linearIndex3D_full) = weightMtx_fullfield(linearIndex3D_full) + (ampCTF(indexingCTF)>0.05);
        end

        if useTimer
            myReconstructionTimer.timeRemain(idx);
        end
        
    end

    mask_fullfield = (weightMtx_fullfield>10^-3);
    recScatteringPotentialFT_fullfield = recScatteringPotentialFT_fullfield./(weightMtx_fullfield + 10^-5);
    
    if useTimer
        myReconstructionTimer.delete;
    end
    
    if valueOnHold
        ftRecons = ftRecons + ftExpanded./(maskExpanded + eps);
        maskRecons = maskRecons + (maskExpanded ~= 0);
    end
    
end

