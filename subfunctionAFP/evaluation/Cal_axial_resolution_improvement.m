function Ratio = Cal_axial_resolution_improvement(NA_illu,NA_illu_extend,n_media)
    % This function calculates the theoretical axial resolution improvement
    % with darkfield extension compared with NA-matching only
    % By Zhenyu Dong
    
    Rz1 = n_media-sqrt(n_media^2-NA_illu^2);
    Rz2 = n_media-sqrt(n_media^2-NA_illu_extend^2);
    Ratio = Rz2/Rz1;
end