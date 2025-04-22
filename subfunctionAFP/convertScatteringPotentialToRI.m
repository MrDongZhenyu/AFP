function [RI_3D_real,RI_3D_imag] = convertScatteringPotentialToRI(scattering_potential, wavelength, n_media)
% Convert scattering potential to refractive index       
% By Zhenyu Dong
      
      k = 2*pi/wavelength;
      alpha = real(scattering_potential)/(k^2);
      beta = imag(scattering_potential)/(k^2);
      RI_3D_real = sqrt((n_media^2+alpha+sqrt((n_media^2+alpha).^2+beta.^2))/2);
      RI_3D_imag = beta/2./RI_3D_real;
end