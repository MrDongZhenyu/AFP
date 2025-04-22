function y = linear_weight(x, x0, x1, offset)
% Calculate weight function with a linear drop
% when x=x0, y=1
% when x=x1, y=offset
% offset controls the slope
% By Zhenyu Dong
    
y = (1-offset)*(x1-x)/(x1-x0) + offset;

end
