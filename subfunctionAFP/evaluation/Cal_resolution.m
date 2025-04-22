function [Rxy,Rz] = Cal_resolution(NA,lambda,n)
    % Calculate the resolution of AFP in 3D 
    % lambda in um
    % By Zhenyu Dong

    Rxy = lambda/2/NA;  % in um
    Rz = lambda/(n-sqrt(n^2-NA^2)); % in um
end