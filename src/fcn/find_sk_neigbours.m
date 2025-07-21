function [sk_neighbours, A] = find_sk_neigbours(sk_structure,...
    branch_points, tip_points)
% Returns list of neighbours for each BP and corresponding morphological
% graph.

n_bps = size(branch_points, 1);
sk_neighbours = {n_bps};
[n, m] = size(sk_structure);
% Create a mask for removing all BPs
mask = zeros(n,m);
idx = sub2ind(size(sk_structure), branch_points(:, 2),...
    branch_points(:, 1));
mask(idx) = 1;
mask = imdilate(mask, strel('square', 3));

sk_structure_no_bps = sk_structure & ~mask;


morph_nodes = [branch_points; tip_points];


n_morph_nodes = size(morph_nodes, 1);

A = zeros(n_morph_nodes);

for i = 1:n_bps
    % Identify neigbors
    current_bp = branch_points(i, :);
    branch_points_tmp = branch_points;
    branch_points_tmp(i, :) = [];
    mask_tmp = zeros(n,m);
    mask_tmp(current_bp(2), current_bp(1)) = 1;
    mask_tmp = imdilate(mask_tmp, strel('square', 3));
    sk_structure_tmp = sk_structure_no_bps | mask_tmp;
    % dist_map = bwdistgeodesic(sk_structure_tmp, current_bp(1), current_bp(2));
    disconnected = imfill(~sk_structure_tmp,...
        [current_bp(2), current_bp(1)], 8);
    sk_structure_tmp = sk_structure_tmp & disconnected;
    sk_structure_tmp = bwskel(sk_structure_tmp); 
    tips = bwmorph(sk_structure_tmp, "endpoints");

    [tip_x, tip_y] = find(tips);
    
    [K, D] = dsearchn(branch_points_tmp, [[tip_y, tip_x]; current_bp]);

    close_bp = (D<=sqrt(18)); % I.e. a diagonal of 3 pixels
    neighbor_coords = branch_points_tmp(K(close_bp), :);
    neighbour_x = neighbor_coords(:, 1);
    neighbour_y = neighbor_coords(:, 2);

    % Compare to tip-points
    [K, D] = dsearchn(tip_points, [tip_y, tip_x]);
    same_tip = (D <=sqrt(18));
    connected_tip_coords = tip_points(K(same_tip), :);
    connected_tip_indices = K(same_tip);
    tip_node_indices = connected_tip_indices + n_bps;
    
    % Calculate geodesic distance to neigbors.
    D = bwdistgeodesic(logical(sk_structure), ...
        current_bp(1), current_bp(2), 'quasi-euclidean');
    idx = sub2ind(size(sk_structure), [neighbour_y;...
        connected_tip_coords(:, 2)], [neighbour_x;...
        connected_tip_coords(:, 1)]);
    bp_distances = D(idx);
    idx = sub2ind(size(sk_structure), connected_tip_coords(:, 2),...
        connected_tip_coords(:, 1));
    tip_distances = D(idx);
    [neighbor_coords, ia, ~] = unique(neighbor_coords, 'rows', 'stable');
    [~, neighbor_indices] = ismember(neighbor_coords, branch_points, 'rows');
    sk_neighbours{i} = {current_bp, neighbor_coords, bp_distances(ia)};
    
    % Morphological graph

    A([neighbor_indices; tip_node_indices] , i) =...
        [bp_distances(ia); tip_distances];
    A(i, [neighbor_indices; tip_node_indices]) = ...
        [bp_distances(ia)' tip_distances'];


end

end