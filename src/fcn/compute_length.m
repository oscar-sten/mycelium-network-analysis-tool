function length = compute_length(sk_structure, mm_per_pixel)
% Length auto should be 4-connected pixels in sk_structure * mm_per_pixel +
% exclusively 8-connected pixels in sk_structure * mm_per_pixel * sqrt(2).

% Input: sk_structure, a skeletonized BW image.
%        mm_per_pixel, the length to which a pixel corresponds.
% Output: length

% All pixels
sum_sk_structure_all = sum(sk_structure(:));

% Identify 4-connected pixels
sk_structure_4_conn = bwareaopen(sk_structure, 2, 4);
sum_sk_structure_4_conn = sum(sk_structure_4_conn(:));

% Compute number of 8-connected pixels
sum_sk_structure_8_conn = sum_sk_structure_all - sum_sk_structure_4_conn;

% Compute length multiplying diagonal pixels with sqrt(2)
length = (sum_sk_structure_4_conn * mm_per_pixel) +...
    (sum_sk_structure_8_conn * mm_per_pixel * sqrt(2));

end