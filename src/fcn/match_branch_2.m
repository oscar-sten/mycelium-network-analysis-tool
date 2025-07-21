function [branch_out_coords, direction_vectors] = match_branch_2(patch,...
    bp, feature_vector, connection_points, x_range, y_range,...
    direction_vectors, theta_swich)

% Identifies the branch out arm in a set of three branches.
% 
% Input: 
%       patch : [X x Y], binary matrix surrounding the BP.
%       feature_vector : [2 x 3], ridge and angle features.
%       connection_points : the points where the arms are connected to the
%                           BP.
%       x_range : [1 x X], x coordinates in patch.
%       y_range : [1 x Y], y coordinates in patch.
%       direction_vectors : [2 x 3], the vectors repressenting the arms'
%                           directions in arbitrary order.
% 
% Output: 
%       branch_out_coords : [1 x 2], the connection coorinates of the
%                           branch out arm.
%       direction_vectors : [2 x 3], the vectors repressenting the arms'
%                           directions. In columns 1 and 2, the main
%                           direction vectors. In column 3, the branch out
%                           vector.



if size(direction_vectors, 2)==3
    point_combinations = nchoosek(1:3, 2);
    phi_list = zeros(1, 3);
    for i=1:3
        pair = point_combinations(i, :);
        vector1 = direction_vectors(:, pair(1));
        vector2 = direction_vectors(:, pair(2));
        dp = dot(vector1, vector2);
        abs_1 = norm(vector1);
        abs_2 = norm(vector2);
        cos_phi = dp/(abs_1*abs_2);
        phi = rad2deg(acos(cos_phi));
        phi_list(i) = phi;
    end
    % [~, phi_max_index] = max(phi_list);
    close_to_right_angles = (theta_swich > abs(phi_list - 90));
    n_close_to_right_angles = sum(close_to_right_angles);
    
    if n_close_to_right_angles == 2
        [~, phi_max_index] = max(phi_list);
        main_indices = point_combinations(phi_max_index, :);
        branch_index = setdiff((1:3)', main_indices');
    else
        [~, branch_index] = min(feature_vector(1, :));
    end

  
else
    [~, branch_index] = min(feature_vector(1, :));
end



branch_out_coords = connection_points(:, branch_index)... 
    + [x_range(1)-1; y_range(1)-1];

manhattan_diff = branch_out_coords - bp';

test = abs(manhattan_diff)>1;
if test(1) || test(2)
    % Find point between branch out coordinates and 
    D_tmp1 = bwdistgeodesic(patch, connection_points(1, branch_index),...
        connection_points(2, branch_index));
    D_tmp2 = bwdistgeodesic(patch, bp(1)-x_range(1)+1, ...
        bp(2)-y_range(1)+1);
    D = D_tmp1 + D_tmp2;
    D = round(D * 8) / 8;
    D(isnan(D)) = inf;
    middle_path = imregionalmin(D);
    middle_path(bp(2)-y_range(1)+1, ...
        bp(1)-x_range(1)+1) = 0;
    middle_path(connection_points(2, branch_index),...
        connection_points(1, branch_index)) = 0;
    [y, x] = find(middle_path);
    if numel(y) == 1
        branch_out_coords = [x + x_range(1)-1; y + y_range(1) - 1];
    else
        path_points = [x, y];
        [K,D] = dsearchn(path_points, [bp(1)-x_range(1)+1, ...
            bp(2)-y_range(1)+1]);

        [~, ind] = min(D);
        ind_closest = K(ind);
   
        point = path_points(ind_closest, :);
        branch_out_coords = point'...
            + [x_range(1)-1; y_range(1)-1];
    end
    branch_out_vector = direction_vectors(:, branch_index);
    direction_vectors(:, branch_index) = [];
    direction_vectors = [direction_vectors branch_out_vector];
end




end