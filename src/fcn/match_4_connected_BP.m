function under_passage_points = match_4_connected_BP(bp, ...
    sk_structure, patch_range)

[n, m] = size(sk_structure);

if size(bp, 1) ==1
    center_point = bp;
else
    center_point = round(mean(bp, 1));
end

% Define patch indices
x_range = center_point(1)-patch_range(1):center_point(1)+patch_range(2);
y_range = center_point(2)-patch_range(1):center_point(2)+patch_range(2);


% Remove invalid indices
x_range(x_range<1)=[];
x_range(x_range>m)=[];
y_range(y_range<1)=[];
y_range(y_range>n)=[];


% Extract patch
patch = sk_structure(y_range, x_range);

% Compute the corresponding coordinates of the branch points inside the
% patch.

in_patch = bp - [x_range(1) y_range(1)] + [1 1];

in_patch_center = mean(bp, 1) - [x_range(1) y_range(1)] + [1 1];

bpImg = bwmorph(patch, 'branchpoints');

disconnected = imfill(~bpImg, [in_patch(:, 2) in_patch(:, 1)], 8);

% Remove disconnected components
bpImg = bpImg & disconnected;

patch2 = patch & ~bpImg;

% Remove pixels between bps
if size(bp, 1) == 2
    test = imdilate(bpImg, strel('square', 3));
    patch2 = patch2 & ~test;
end


if sum(bpImg(:))>1
    [bp_coords_X, bp_coords_Y]  = find(bpImg);
    in_patch = [mean(bp_coords_X), mean(bp_coords_Y)];
end

if size(bp, 1) == 1
    [X_coords_patch, Y_coords_patch] = find(patch2);

    coordinates = [X_coords_patch, Y_coords_patch];

    % Calculate the Euclidean distance between the seed and each coordinate
    distances = pdist2(in_patch, coordinates, 'euclidean');

    % Sort the distances in ascending order and get the indices
    [~, indices] = sort(distances);

    % Get the four closest coordinates
    closestCoordinates = coordinates(indices(1:4), :);
else
    [L,n] = bwlabel(patch2, 8);
    if n~=4
        error(strcat("Wrong number of components. Expected 4, found: ", ...
            string(n)))
    end
    coordinates = zeros(n, 2);
    for i=1:n
        [X_coords_patch, Y_coords_patch] = find(L==i);
        coords_tmp = [X_coords_patch, Y_coords_patch];
        distances = pdist2(in_patch, coords_tmp, 'euclidean');
        [~, indices] = sort(distances);
        closestCoordinates = coords_tmp(indices(1), :);
        coordinates(i, :) = closestCoordinates;
    end
    closestCoordinates = coordinates;
end



% Select one pair of coordinates, then match it with the point with the
% furthest geodesic Manhattan distance in a frame where the BP(s) is/are
% excluded
if size(bp, 1) == 1
    notBPPoints = ones(size(patch))&~bpImg;
else
    notBPPoints = ones(size(patch))&~test;
end

seedPoint = closestCoordinates(1, :);

D = bwdistgeodesic(notBPPoints, seedPoint(2), seedPoint(1), 'cityblock');

dists = zeros(3, 1);
for i = 1:3
    dists(i) = D(closestCoordinates(i+1, 1), closestCoordinates(i+1, 2));
end

% Find the index that maximizes the geodesic Manhattan 
[~, max_dist_ind] = max(dists);


under_passage_points = [closestCoordinates([1, max_dist_ind+1], 2), ...
    closestCoordinates([1, max_dist_ind+1], 1)]+...
    [x_range(1) y_range(1)] - [1 1];


end