function [candidate_pairs, very_close_pairs, removed_pairs]...
    = find_possible_crossing_pairs(...
    G, direction_vector_sets_3, branch_points, theta_one)

candidate_pairs = [];
very_close_pairs = [];
removed_pairs = [];
n_vector_sets = size(direction_vector_sets_3, 3);

for i=1:n_vector_sets
    connected_BPs = neighbors(G, i);
    connected_BPs(connected_BPs>n_vector_sets) = [];
    connected_BP_distances = distances(G, i, connected_BPs);

    if  ~isempty(connected_BPs) && (isempty(candidate_pairs)...
            || ~ismember(i, candidate_pairs(:, 2)))
        
        [closest_dist, closest_dist_ind] = min(connected_BP_distances);
        cand_neigbour = connected_BPs(closest_dist_ind);
        neighbours_neigbor = neighbors(G, cand_neigbour);
        neighbours_neighbour_dist = distances(G, cand_neigbour,...
            neighbours_neigbor);
        [~, neighbours_neigbor_closest_neighbor_index] =...
            min(neighbours_neighbour_dist);

        if closest_dist <= 2*sqrt(2)
                 % If BPs have a distance <= 2 diagonals 
                 very_close_pairs = [very_close_pairs; [i cand_neigbour]];
        elseif  neighbours_neigbor(neighbours_neigbor_closest_neighbor_index)...
            == i

            phi = rad2deg(compute_angle(direction_vector_sets_3(:,3, i), ...
                direction_vector_sets_3(:,3, cand_neigbour)));

            A = [1 0;
                0 -1]; % Transfromation matrix

            % Remove points whose connection is in the same direction as
            % the vector.
            between_bp_vector_1 = branch_points(cand_neigbour, :) - ...
                branch_points(i, :);
            between_bp_vector_1 = A * between_bp_vector_1';
            % First way
            phi_2 = rad2deg(compute_angle(between_bp_vector_1,...
                direction_vector_sets_3(:,3, i)));
            % Second way
            between_bp_vector_1 = branch_points(i, :) -...
                branch_points(cand_neigbour, :);
            between_bp_vector_1 = A * between_bp_vector_1';
            phi_3 = rad2deg(compute_angle(between_bp_vector_1,...
                branch_points(cand_neigbour, :)));

       
            if phi > theta_one && (phi_2 > theta_one && phi_3 > theta_one)
                candidate_pairs = [candidate_pairs; [i cand_neigbour]];
            else
                removed_pairs = [removed_pairs; [i cand_neigbour]];
            end
        end
    end
end

end