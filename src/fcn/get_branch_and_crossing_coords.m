function [branch_points, connection_points, under_passage_points,...
    crossing_centers] =...
    get_branch_and_crossing_coords(branch_points, connection_points,...
    perfect_crossings_coords,...
    crossing_pairs, perfect_crossing_under_passage_coords,...
    very_close_pairs, very_close_crossing_under_passage_coords)


%% Structure crossings

if ~isempty(very_close_pairs)
    crossing_centers_x = [mean([branch_points(crossing_pairs(:,1), 1),...
        branch_points(crossing_pairs(:,2), 1)], 2);...
        mean([branch_points(very_close_pairs(:,1), 1),...
        branch_points(very_close_pairs(:,2), 1)], 2);...
        perfect_crossings_coords(:, 1)];

    crossing_centers_y = [mean([branch_points(crossing_pairs(:,1), 2),...
        branch_points(crossing_pairs(:,2), 2)], 2);...
        mean([branch_points(very_close_pairs(:,1), 2),...
        branch_points(very_close_pairs(:,2), 2)], 2);...
        perfect_crossings_coords(:, 2)];
else
    crossing_centers_x = [mean([branch_points(crossing_pairs(:,1), 1),...
        branch_points(crossing_pairs(:,2), 1)], 2);...
        perfect_crossings_coords(:, 1)];

    crossing_centers_y = [mean([branch_points(crossing_pairs(:,1), 2),...
        branch_points(crossing_pairs(:,2), 2)], 2);...
        perfect_crossings_coords(:, 2)];
end

crossing_centers = [crossing_centers_x, crossing_centers_y];

n_matched_crossings = size(crossing_pairs, 1);
n_perfect_crossings = size(perfect_crossings_coords, 1);
n_very_close_crossings = size(very_close_pairs, 1);


matched_crossings = [crossing_centers(1:n_matched_crossings, :)...
    connection_points(crossing_pairs(:, 1), :)...
    connection_points(crossing_pairs(:, 2), :)];

perfect_crossings = [perfect_crossings_coords ...
    perfect_crossing_under_passage_coords(1:2:2*n_perfect_crossings-1, :)...
    perfect_crossing_under_passage_coords(2:2:2*n_perfect_crossings, :)];

if ~isempty(very_close_pairs)
    very_close_crossing_centers = [mean([branch_points(...
        very_close_pairs(:, 1),  1), ...
        branch_points(very_close_pairs(:, 2),  1)], 2),...
        mean([branch_points(very_close_pairs(:, 1),  2), ...
        branch_points(very_close_pairs(:, 2),  2)], 2)];

    very_close_crossings = [very_close_crossing_centers,...
        very_close_crossing_under_passage_coords(1:2:2*n_very_close_crossings-1, :),...
        very_close_crossing_under_passage_coords(2:2:2*n_very_close_crossings, :)];
else
    very_close_crossings = [];
end


% Merge double "perfect crossings"
double_indices = [];
for i=1:n_perfect_crossings-1
    dist = norm(perfect_crossings_coords(i, :)...
        - perfect_crossings_coords(i+1, :));
    if dist<sqrt(3)
        double_indices = [double_indices; i; i+1];
    end
end

remove = [];
for i=1:numel(double_indices)/2
    pair_ind = double_indices(2*i-1:2*i, :);
    perfect_crossings(pair_ind(1), 1:2) = mean(...
        [perfect_crossings(pair_ind(1), 1:2);...
        perfect_crossings(pair_ind(2), 1:2)], 1);
    
    % Identify the connection points closest to the crossing center and set
    % them as connection points.
    dist1 = norm(perfect_crossings(pair_ind(1), 1:2)...
        -perfect_crossings(pair_ind(1), 3:4));
    dist2 = norm(perfect_crossings(pair_ind(1), 1:2)...
        -perfect_crossings(pair_ind(2), 3:4));
    dist3 = norm(perfect_crossings(pair_ind(1), 1:2)...
        -perfect_crossings(pair_ind(1), 5:6));
    dist4 = norm(perfect_crossings(pair_ind(1), 1:2)... 
        -perfect_crossings(pair_ind(2), 5:6));

   [~, k] = min([dist1, dist2]);

    perfect_crossings(pair_ind(1), 3:4) = perfect_crossings(...
        pair_ind(k), 3:4);

   [~, k] = min([dist3, dist4]);

    perfect_crossings(pair_ind(1), 5:6) = perfect_crossings(...
        pair_ind(k), 5:6);


    remove = [remove; pair_ind(2)];
end

perfect_crossings(remove, :) = [];

crossings = [matched_crossings; perfect_crossings; very_close_crossings];


under_passage_points = [];
for i=1:size(crossings, 1)
    under_passage_points = [under_passage_points; crossings(i, 3:4);...
        crossings(i, 5:6)];
end

%% Remove invalid branch and connection points

branch_points(ismember(branch_points, perfect_crossings_coords,...
    'rows'), :) = [];


% Remove paired crossings from set of branch points
if ~isempty(very_close_pairs)
    non_paired_bp_indices = setdiff(1:size(branch_points, 1), [crossing_pairs(:,1);...
        crossing_pairs(:, 2); very_close_pairs(:, 1); very_close_pairs(:, 2)]);
else
    non_paired_bp_indices = setdiff(1:size(branch_points, 1), [crossing_pairs(:,1);...
        crossing_pairs(:, 2)]);
end

branch_points = branch_points(non_paired_bp_indices, :);
connection_points = connection_points(non_paired_bp_indices, :);


end