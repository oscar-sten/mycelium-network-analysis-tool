function [branch_points, perfect_crossings_coords,  direction_vector_sets_3,...
    direction_vector_sets_4, connection_points] = get_BP_vectors(...
    sk_structure, ridgeFilt, patch_range, theta_swich)


%% Indentify all branch points
[branch_points_Y, branch_points_X] = find(...
    bwmorph(sk_structure,'branchpoints'));

% %% Calculate adjacencies between branch points
branch_points = [branch_points_X branch_points_Y];
[brp_n, brp_m] = size(branch_points);


%% Analyse all BPs as branch points

connection_points = zeros(brp_n, brp_m);
% direction_vector_sets = zeros(2,3,brp_n);
direction_vector_sets_3 = [];
direction_vector_sets_4 = [];
false_bp_indices = [];
perfect_crossings_indices = [];

for i=1:brp_n
    current_bp = branch_points(i, :);
    [connection_point, ~, ~, direction_vectors]...
    = single_bp_analysis(current_bp, sk_structure, ridgeFilt,...
    patch_range, theta_swich);
   connection_points(i, :) = connection_point;
   if numel(direction_vectors)<6
        false_bp_indices = [false_bp_indices; i];
   elseif numel(direction_vectors) > 6
        perfect_crossings_indices = [perfect_crossings_indices; i];
        direction_vector_sets_4 = cat(3, ...
            direction_vector_sets_4, direction_vectors);
   else
       direction_vector_sets_3 = cat(3, ...
           direction_vector_sets_3, direction_vectors);
   end

end


%% Remove unwanted indices

perfect_crossings_coords = branch_points(perfect_crossings_indices, :);
branch_points([false_bp_indices; perfect_crossings_indices], :) = [];

connection_points([false_bp_indices; perfect_crossings_indices], :) = [];



branch_points = [branch_points; perfect_crossings_coords];


end