function [crossing_pairs, removed_pairs_step_2] = ...
    find_likely_crossing_pairs(G, direction_vector_sets_3,...
    candidate_pairs, branch_points, dist_limit, theta_T)

removed_pairs_step_2 = [];
crossing_pairs = [];

n_candidate_pairs = size(candidate_pairs, 1);

for i=1:n_candidate_pairs
    pair_indices = candidate_pairs(i, :);
    pair = branch_points(pair_indices, :);
    dist = distances(G, pair_indices(1), pair_indices(2));
    % Compute branch vector angle of first bp
    bp_vector_set_1 = direction_vector_sets_3(:, :, pair_indices(1));
    main_vector = bp_vector_set_1(:, 1) - bp_vector_set_1(:, 2);
    branch_vector = bp_vector_set_1(:, 3);
    dp = dot(main_vector, branch_vector);
    abs_1 = norm(main_vector);
    abs_2 = norm(branch_vector);
    cos_phi = dp/(abs_1*abs_2);
    phi_bp_1 = rad2deg(acos(cos_phi));
    % Compute branch vector angle of second bp
    bp_vector_set_2 = direction_vector_sets_3(:, :, pair_indices(2));
    main_vector = bp_vector_set_2(:, 1) - bp_vector_set_2(:, 2);
    branch_vector = bp_vector_set_2(:, 3);
    dp = dot(main_vector, branch_vector);
    abs_1 = norm(main_vector);
    abs_2 = norm(branch_vector);
    cos_phi = dp/(abs_1*abs_2);
    phi_bp_2 = rad2deg(acos(cos_phi));
    

    is_crossing = crossing_criterion(dist, dist_limit, phi_bp_1,...
        phi_bp_2, theta_T);
    if is_crossing
        % Remove from candidate set and add to removed
        crossing_pairs = [crossing_pairs; candidate_pairs(i, :)];
    else
        removed_pairs_step_2 = [removed_pairs_step_2;...
            candidate_pairs(i, :)];
    end
end


end