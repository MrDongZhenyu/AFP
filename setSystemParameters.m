function p = setSystemParameters(sampleName)
% This function sets the parameters of the imaging system and the sample
% By Zhenyu Dong

    p.mag = 20;                   % magnification of the objective lens
    if sampleName == "bead"
        % LED array setup
        p.ps = 3.45/p.mag;        % unit: micron. Pixel size in object space
        p.NA = 0.412;             % numerical aperture (NA) of the imaging system
        p.lambda = 520;           % unit: nm. wavelength of the illumination light
        p.useAbsorption = false;  % whether to generate absorption image

        p.n_media = 1.58;         % refractive index of the environment (RI liquid)
        p.pixelSizeZ = 2*p.ps;    % unit: micron. pixelSizeZ can be 
                                  % different from pixelSizeXY
        p.Diameter = 6;           % unit: micron.
        p.samThick = min(round(p.Diameter/p.pixelSizeZ*1.5),round(50/p.pixelSizeZ)); % unit: pixel.
        p.zsize    = 100;         % unit: pixel.
        p.nDFperR = 48*ones(1,3); % Darkfield illumination number per ring.

    elseif sampleName == "embryo"
        % Laser setup
        p.ps = 6.5/p.mag;         % unit: micron. Pixel size in object space
        p.NA = 0.40;              % numerical aperture (NA) of the imaging system
        p.lambda = 532;           % unit: nm. wavelength of the illumination light
        p.useAbsorption = false;  % whether to generate absorption image

        p.n_media = 1.33;         % refractive index of the environment (KSOM or PBS)
        p.pixelSizeZ = 1;         % unit: micron. pixelSizeZ can be 
                                  % different from pixelSizeXY
        p.zsize  = 105;           % unit: pixel.
        p.samThick = p.zsize*0.4; % unit: pixel.
        p.nDFperR = [45,48,52,55,58,60]; % Darkfield illumination number per ring.
    
    elseif sampleName == "root"
        % LED array setup
        p.ps = 3.45/p.mag;        % unit: micron. Pixel size in object space
        p.NA = 0.412;             % numerical aperture (NA) of the imaging system
        p.lambda = 520;           % unit: nm. wavelength of the illumination light
        p.useAbsorption = false;  % whether to generate absorption image

        p.n_media = 1.40;         % refractive index of the environment (glycerol & water, 1:1)
        p.pixelSizeZ = 0.7;       % unit: micron. pixelSizeZ can be 
                                  % different from pixelSizeXY
        p.zsize  = 200;           % unit: pixel.
        p.samThick = p.zsize*0.35;% unit: pixel.
        p.nDFperR = 48*ones(1,6); % Darkfield illumination number per ring.

    elseif sampleName == "pathology"
        % DUV setup
        p.ps = 2.74/p.mag;        % unit: micron. Pixel size in object space
        p.NA = 0.36;              % numerical aperture (NA) of the imaging system
        p.lambda = 265;           % unit: nm. wavelength of the illumination light
        p.useDarkfield = false;   % whether to use darkfield reconstruction
        p.useAbsorption = true;   % whether to generate absorption image

        p.n_media = 1.4998;       % refractive index of the environment (crystal mount)
        p.pixelSizeZ = 0.8;       % unit: micron. pixelSizeZ can be 
                                  % different from pixelSizeXY
        p.zsize  = 60;            % unit: pixel.
        p.samThick = p.zsize*0.4; % unit: pixel.
    
    else
        error('Please define your imaging system and sample parameter first!'); 
        % Add you imaging system and sample parameters here.
    end

    p.pixelSizeXY = p.ps;  % unit: micron. This is the effective pixel size (lateral)
  
end