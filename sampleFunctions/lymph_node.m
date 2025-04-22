function obj = lymph_node(target_size,RI_bg,SamThick)
    % Generate lymph node phantom
    % By Haowen Zhou, Ruizhi Cao

    if nargin < 3 || isempty(SamThick)
        SamThick = 50; 
    end
    if nargin < 2 || isempty(RI_bg)
        RI_bg = 1.33068; 
    end
    if nargin < 2 || isempty(target_size)
        target_size = [200,200,200]; 
    end

    if SamThick > target_size(3)
        disp('Warning zsize is smaller than number of sample slices (SamThick). SamThick is set to round(zsize/2)');
        SamThick = round(target_size(3)/2);
    end

    load('lymphn_cropped_8bit.mat');
    [xsizeObj,ysizeObj,zsizeObj] = size(obj);

    % If the sample is a non-square matrix in lateral dimensions 
    if xsizeObj ~= ysizeObj
        size2use = max(xsizeObj,ysizeObj);
        padsize = [size2useObj - xsizeObj,size2use - ysizeObj];
        obj = padarray(obj,padsize,RI_bg,'post');
        warning('This code assumes the lateral dimensions of the object are the same. The object will be zero padded to satisfy this constraint.');
    end

    % Pad to the target size for xy dimensions
    [xsizeObj,ysizeObj,~] = size(obj);
    if xsizeObj < target_size(1)
        padsize = [floor((target_size(1) - xsizeObj)/2),floor((target_size(2) - ysizeObj)/2), 0];
        if mod(target_size(1) - xsizeObj,2) == 0
            obj = padarray(obj,padsize,RI_bg,'both');
        else 
            obj = padarray(obj,padsize,RI_bg,'both');
            obj = padarray(obj,[1,1,0],RI_bg,'post');
        end
    end

    % Select target number of z-slices
    if zsizeObj > SamThick
        centerZ = ceil(zsizeObj/2);
        obj = obj(:,:,centerZ-floor(SamThick/2):centerZ+floor(SamThick/2)-1);
    end
    
    % Pad to the target number of z-slices
    [~,~,zsizeObj] = size(obj);
    if zsizeObj < target_size(3)
        padsize = floor((target_size(3) - zsizeObj)/2);
        if mod(target_size(3) - zsizeObj,2) == 0
            obj = padarray(obj,[0,0,padsize],RI_bg,'both');
        else 
            obj = padarray(obj,[0,0,padsize],RI_bg,'both');
            obj = padarray(obj,[0,0,1],RI_bg,'post');
        end
    end

%     % Check sample volume
%     figure;
%     for idx = 1:target_size(3)
%         imagesc(obj(:,:,idx));axis square;clim([RI_bg,RI_bg+0.01]);title(num2str(idx));
%         pause (0.1);
%     end

end