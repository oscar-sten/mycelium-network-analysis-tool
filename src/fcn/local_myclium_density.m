function [local_density_map, global_density] = local_myclium_density(...
    sk_structure, mm_per_pixel)
% Returns local mycelium density in mm/mm^2

[n, m] = size(sk_structure);
global_density =  sum(sk_structure(:))/(n*m*mm_per_pixel);
% window = round(368000 * mm_per_pixel); % Design parameter, 400 pixels
window = round(460000 * mm_per_pixel); % Design parameter, 500 pixels
M = filter2(ones(window), sk_structure);
local_density_map = M / (window * window * mm_per_pixel);
end