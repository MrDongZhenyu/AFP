function DOF = Cal_DOF(NA,lambda,n,Mag,ps)
    % This function calculates the Depth of field (DOF) of a microscope
    % lambda and ps in um
    % n is the refractive index of the mounting media
    % NA: numerical aperture of the objective lens
    % Mag: magnification of the objective lens
    % ps is the pixel size of the camera sensor
    % By Zhenyu Dong

    DOF = lambda*n/(NA)^2 + n*ps/Mag/NA; % in um
end