function [isolated_filaments, hyphal_labels, faulty_crossings] = ...
    isolate_filaments(sk_structure,...
    all_branches, under_cross_coords, sorted)

%% Separates the filaments given a binary skeleton 


sk_st_isolated_hyphae = sk_structure;
img_dimension = size(sk_structure);

xt = all_branches(:, 1);
yt = all_branches(:, 2);
idx = sub2ind(img_dimension, yt, xt);


%% Remove branches from binary skeleton
sk_st_isolated_hyphae(idx) = 0;


%% Remove underlaying crossing points from binary skeleton
xt = under_cross_coords(:, 1);
yt = under_cross_coords(:, 2);
idx = sub2ind(img_dimension, yt, xt);
sk_st_isolated_hyphae(idx) = 0;

[lbl_sk_st_isolated_hyphae, n_isolated_hyphae] = bwlabel(...
    sk_st_isolated_hyphae, 8);

%% Sort under_cross coordinates.
% Assert that number of points is even.
[n_under_coords, ~] = size(under_cross_coords);
assert(rem(n_under_coords, 2) == 0);

% Build adjacency matrix with Euclidean or Geodesic distances between all points.
under_cross_adj_mat = zeros(n_under_coords);
for i=1:n_under_coords
    point_coords = under_cross_coords(i, :);
    for j=1:n_under_coords
        if i~=j
            dist = norm(point_coords-under_cross_coords(j, :));
        else
            dist = Inf;
        end
        under_cross_adj_mat(i,j) = dist;
        under_cross_adj_mat(j, i) = dist;
    end
end

% While non Inf elements exist in adjacency matrix: 
% Remove identify the smallest value in the adjacency matrix.
% Add the identified points to the list.
% Set the respective rows/columns to Inf
if ~sorted
    under_cross_coords_sorted = [];
    while sum(isinf(under_cross_adj_mat(:))) < numel(under_cross_adj_mat)
        [~,I] = min(under_cross_adj_mat,[],'all');
        [row,col] = ind2sub([n_under_coords n_under_coords], I);
        under_cross_adj_mat(row, :) = Inf;
        under_cross_adj_mat(col, :) = Inf;
        under_cross_adj_mat(:, row) = Inf;
        under_cross_adj_mat(:, col) = Inf;
        under_cross_coords_sorted = [under_cross_coords_sorted;
            under_cross_coords(row, :);
            under_cross_coords(col, :)];
    end
else
    under_cross_coords_sorted = under_cross_coords;
end


%% Identfy IDs of filaments connected to eachother by undergoing coordinates.

% Match under-crossing coordinates
nodeIDCrossingPairs = [];
faulty_crossings = [];
% Loop over pairs of under-going coordinates
for i=1:n_under_coords/2
    pair = under_cross_coords_sorted(2*i-1:2*i, :);
    % Extract the surrounding coordinates of both points. Since both of
    % them border with the "master hypha" the label of that hypha will be
    % shared between them. They will have one unique number label each
    % which is the number of their associated underlying hypha.
    test_patch_1 = lbl_sk_st_isolated_hyphae(...
        pair(1, 2)-1:pair(1, 2)+1,...
        pair(1, 1)-1:pair(1, 1)+1);
    test_patch_2 = lbl_sk_st_isolated_hyphae(...
        pair(2, 2)-1:pair(2, 2)+1,...
        pair(2, 1)-1:pair(2, 1)+1);
    test_1 = nonzeros(test_patch_1);
    test_2 = nonzeros(test_patch_2);
    first_unique = setdiff(test_1,test_2);
    second_unique = setdiff(test_2,test_1);
    % if ~isempty(first_unique) && ~isempty(second_unique)
    if (numel(first_unique) == 1) && (numel(second_unique) == 1)
        nodeIDCrossingPairs = [nodeIDCrossingPairs;...
            [first_unique second_unique]];
        % If this operation is not possible, it means that there is not a
        % unique matching available and an error has occured earlier.
    elseif numel(first_unique) > 1 || numel(second_unique) > 1
        faulty_crossings = [faulty_crossings; pair];
    elseif isempty(first_unique) && isempty(second_unique)
        nodeIDCrossingPairs = [nodeIDCrossingPairs;...
            [test_1(1) test_1(2)]];
    elseif isempty(first_unique)
        nodeIDCrossingPairs = [nodeIDCrossingPairs;...
            [test_1(1) second_unique]];
    elseif isempty(second_unique)
        nodeIDCrossingPairs = [nodeIDCrossingPairs;...
            [first_unique test_2(1)]];
    end
end

% Sort rows to have always right order
nodeIDCrossingPairs = sort(nodeIDCrossingPairs, 2);
% The below row may be symptom of a problem
nodeIDCrossingPairs = unique(nodeIDCrossingPairs, 'rows');
back_tracked = 0;
hyphal_labels = 1:n_isolated_hyphae;


%% Loop over the set of filament IDs matched to each other and set
% connected filaments to have equal IDs. Algorithm 1.

traversals = traversalsFromPairs(nodeIDCrossingPairs);

for i=1:numel(traversals)
    traversal_ids = traversals{i};
    origin_id = traversal_ids(1);
    connected_ids = traversal_ids(2:end);
    for j=1:numel(connected_ids)
        lbl_sk_st_isolated_hyphae(lbl_sk_st_isolated_hyphae==...
            connected_ids(j)) = origin_id;
    end
    hyphal_labels(connected_ids) = 0;
end


hyphal_labels(hyphal_labels==0) = [];
isolated_filaments = lbl_sk_st_isolated_hyphae;

end