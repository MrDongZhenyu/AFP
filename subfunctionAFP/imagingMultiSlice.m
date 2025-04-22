function [imStack] = imagingMultiSlice(obj,nNAmatching,k_illu,CTF,wavelength,pixelXY,pixelZ,varargin)
% generate images of an extended 3D object using the multi-slice beam
% propagation method
% By Ruizhi cao, Zhenyu Dong

outputField = false;   % whether to output complex field or intensity
useImsize   = false;   % whether to use user defined image size
n_media     = 1.33068; % refractive index of the background media
pad_size    = 0;       % no padding in default
use_gpu     = false;   % no GPU in default

if ~isempty(varargin)
    idx = 1;
    while idx <= length(varargin)
        switch lower(varargin{idx})
            case 'field'
                outputField = true;
                idx = idx+1;
            case 'intensity'
                outputField = true;
                idx = idx+1;
            case {'imsize','image size'}
                imsize = varargin{idx+1};
                useImsize = true;
                idx = idx+2;
            case {'n','refractive index'}
                n_media = varargin{idx+1};
                idx = idx+2;
            case {'center slice','centralslice','central slice','slice center'}
                numSliceToCenter = varargin{idx+1};
                idx = idx+2;
            case {'pad','pad_size'}
                pad_size = varargin{idx+1};
                idx = idx+2;
            case {'use_darkfield','darkfield'}
                use_darkfield = varargin{idx+1};
                idx = idx+2;
            case {'use_gpu','gpu'}
                use_gpu = varargin{idx+1};
                idx = idx+2;
            otherwise
                error('Unsupported option.');
        end
    end
end

obj = padarray(obj,[pad_size,pad_size,0],n_media);  % padding object in x-y domain
lambda = wavelength*10^(-3);                        % unit: micron
[xsizeObj,ysizeObj,zsizeObj] = size(obj);
if useImsize
    xsize = imsize(1);
    ysize = imsize(end);
else
    xsize = xsizeObj;
    ysize = ysizeObj;
end
pixelDownsamp = xsizeObj/xsize;
nIllumination = numel(k_illu)/2;

[YObj,XObj] = meshgrid(1:ysizeObj,1:xsizeObj);
xcObj = floor(xsizeObj/2+1);
ycObj = floor(ysizeObj/2+1);
RObj = abs((XObj-xcObj) + 1i*(YObj-ycObj));
bdCrop = calBoundary([xcObj,ycObj],[xsize,ysize]);

kxy = RObj/pixelDownsamp/(xsize*pixelXY);

% construct "free space propagator" for each slice
zProp_FT = exp(1i*2*pi*sqrt((n_media/lambda)^2-kxy.^2)*pixelZ);

if ~exist('numSliceToCenter','var')
    numSliceToCenter = floor(zsizeObj/2);
end

% allocate space for the results
if use_darkfield
    illu_num = nIllumination;
else
    illu_num = nNAmatching;
end
imStack = zeros(xsize-2*pad_size,ysize-2*pad_size,illu_num);

if use_gpu
    obj = gpuArray(obj);
    zProp_FT = gpuArray(zProp_FT);
    imStack = gpuArray(imStack);
end

for idx = 1:illu_num
    if mod(idx,20)==0
        disp(['Angle: ',num2str(idx)]);
    end
    temp = zeros(xsizeObj,ysizeObj);
    temp(xcObj + round(k_illu(idx,1)*xsize/(xsize-2*pad_size)),ycObj + round(k_illu(idx,2)*ysize/(ysize-2*pad_size))) = xsizeObj*ysizeObj;
    if use_gpu
        temp = gpuArray(temp);
    end
    incField = ifft2(ifftshift(temp));
    
    % forward propagation
    for idxZ = 1:zsizeObj
        outField = incField.*exp(1i*2*pi*( (obj(:,:,idxZ) - n_media)*pixelZ/lambda ));
        outFieldFT = fftshift(fft2(outField));
        incField = ifft2(ifftshift(outFieldFT.*zProp_FT));
    end
    
    % propagte back to focal plane
    field_MS_FT = outFieldFT.*exp(-1i*2*pi*sqrt((n_media/lambda)^2-kxy.^2)*pixelZ*(zsizeObj-1-numSliceToCenter));
    field_MS = ifft2(ifftshift(field_MS_FT(bdCrop(1):bdCrop(2),bdCrop(3):bdCrop(4)).*CTF));
    
    % crop out the paddings
    if outputField
        imStack(:,:,idx) = field_MS((pad_size+1):(end-pad_size),(pad_size+1):(end-pad_size)); 
    else
        temp = field_MS.*conj(field_MS);
        imStack(:,:,idx) = temp((pad_size+1):(end-pad_size),(pad_size+1):(end-pad_size));     
    end
end

if use_gpu
    imStack = gather(imStack);
    imStack = double(imStack);
end

end

