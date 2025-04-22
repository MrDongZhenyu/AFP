function imStack = imagingMultiSliceStack(obj,k_illu,CTF,wavelength,pixelXY,pixelZ,varargin)
% Generate image Stack of an extended object using the multi-slice model
% Input one illumination angle, and output a 3D image Stack, simulating the
% images obtained when we do through-focus scanning through the z-plane.
% Note: k_illu is in pixels, and can be a decimal
% By Zhenyu Dong

outputField = false;   % whether to output complex field or intensity
useImsize   = false;   % whether to use user defined image size
n_media     = 1.33068; % refractive index of the surrounding media
pad_size    = 0;       % no padding in default
use_gpu     = false;   % donot use GPU in default

if ~isempty(varargin)
    idx = 1;
    while idx <= length(varargin)
        switch lower(varargin{idx})
            case 'field'
                outputField = true;
                idx = idx+1;
            case {'imsize','image size'}
                imsize = varargin{idx+1};
                useImsize = true;
                idx = idx+2;
            case {'n','refractive index'}
                n_media = varargin{idx+1};
                idx = idx+2;
            case {'pad','pad_size'}
                pad_size = varargin{idx+1};
                idx = idx+2;
            case {'use_gpu','gpu'}
                use_gpu = varargin{idx+1};
                idx = idx+2;
            otherwise
                error('Unsupported option.');
        end
    end
end

obj = padarray(obj,[pad_size,pad_size,0],n_media); % padding object in x-y domain
lambda = wavelength*10^(-3); % unit: micron
[xsizeObj,ysizeObj,zsizeObj] = size(obj);
if useImsize
    xsize = imsize(1);
    ysize = imsize(end);
else
    xsize = xsizeObj;
    ysize = ysizeObj;
end
pixelDownsamp = xsizeObj/xsize;

[YObj,XObj] = meshgrid(1:ysizeObj,1:xsizeObj);
xcObj = floor(xsizeObj/2+1);
ycObj = floor(ysizeObj/2+1);
RObj = abs((XObj-xcObj) + 1i*(YObj-ycObj));
bdCrop = calBoundary([xcObj,ycObj],[xsize,ysize]);

kxy = RObj/pixelDownsamp/(xsize*pixelXY);

% construct "free space propagator" for each slice
zProp_FT = exp(1i*2*pi*sqrt((n_media/lambda)^2-kxy.^2)*pixelZ);

% calculation start:
imStack = zeros(xsize-2*pad_size,ysize-2*pad_size,zsizeObj);

if use_gpu
    obj = gpuArray(obj);
    zProp_FT = gpuArray(zProp_FT);
    imStack = gpuArray(imStack);
end

% a planewave incidence
temp = zeros(xsizeObj,ysizeObj);
temp(xcObj + round(k_illu(1)*xsize/(xsize-2*pad_size)),ycObj + round(k_illu(2)*ysize/(ysize-2*pad_size))) = xsizeObj*ysizeObj;
if use_gpu
    temp = gpuArray(temp);
end
incField = ifft2(ifftshift(temp));
    
% forward propagation
for idxZ = 1:zsizeObj
    outField = incField.*exp(1i*2*pi*((obj(:,:,idxZ) - n_media)*pixelZ/lambda ));
    outFieldFT = fftshift(fft2(outField));
    incField = ifft2(ifftshift(outFieldFT.*zProp_FT));
end

% propagate back to the image volume
for idxZ = 1:zsizeObj
    field_MS_FT = outFieldFT.*exp(-1i*2*pi*sqrt((n_media/lambda)^2-kxy.^2)*pixelZ*(zsizeObj-idxZ));
    field_MS = ifft2(ifftshift(field_MS_FT(bdCrop(1):bdCrop(2),bdCrop(3):bdCrop(4)).*CTF));
    
    % Crop out the paddings
    if outputField
        imStack(:,:,idxZ) = field_MS((pad_size+1):(end-pad_size),(pad_size+1):(end-pad_size));
    else
        temp = field_MS.*conj(field_MS);
        imStack(:,:,idxZ) = temp((pad_size+1):(end-pad_size),(pad_size+1):(end-pad_size));
    end
end

if use_gpu
    imStack = gather(imStack);
    imStack = double(imStack);
end

end


