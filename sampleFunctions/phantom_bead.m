function phantom = phantom_bead(outer_size,inner_size,RI_bg,RI_bead)
    % Generate bead phantom
    % By Zhenyu Dong

    phantom = double(RI_bg*ones(outer_size));
    d1 = reshape(1:outer_size(1),[],1,1)-(floor(outer_size(1)/2)+1);
    d2 = reshape(1:outer_size(2),1,[],1)-(floor(outer_size(2)/2)+1);
    d3 = reshape(1:outer_size(3),1,1,[])-(floor(outer_size(3)/2)+1);
    
    d1_norm = 2.*d1./inner_size(1);
    d2_norm = 2.*d2./inner_size(2);
    d3_norm = 2.*d3./inner_size(3);
    r_norm = sqrt(d1_norm.^2+d2_norm.^2+d3_norm.^2);
    phantom = phantom + double((RI_bead-RI_bg)*(r_norm<1));
    
end