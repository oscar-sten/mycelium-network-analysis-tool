function G = gen_graph(sk_structure, lbl_sk_st_isolated_hyphae, ...
    all_branches, under_cross_coords, hyphal_labels, ...
    hyphal_median_ridge_values)

% Remove remove objects witout any branches
sk_st_isolated_hyphae = sk_structure;
img_dimension = size(sk_structure);

xt = all_branches(:, 1);
yt = all_branches(:, 2);
idx = sub2ind(img_dimension, yt, xt);

% Remove branches
sk_st_isolated_hyphae(idx) = 0;


% Remove underlaying crossing points 
xt = under_cross_coords(:, 1);
yt = under_cross_coords(:, 2);
idx = sub2ind(img_dimension, yt, xt);
sk_st_isolated_hyphae(idx) = 0;

remaining_isolated_objects = sk_st_isolated_hyphae;

[~, main_hypha_index] = max(hyphal_median_ridge_values);

hypha_id = hyphal_labels(main_hypha_index);

img_dimension = size(sk_structure);
branch_mat = zeros(img_dimension);
xt = all_branches(:, 1);
yt = all_branches(:, 2);
idx = sub2ind(img_dimension, yt, xt);
branch_mat(idx) = 1;

hypha_id_queaue = hypha_id;

Adjacency_matrix = zeros(numel(hyphal_labels), numel(hyphal_labels));
node_names = string(hyphal_labels);


while ~isempty(find(remaining_isolated_objects, 1))
    % Breadth first search
    try
        hypha_id = hypha_id_queaue(1);
    catch
        hypha_id_queaue = unique(lbl_sk_st_isolated_hyphae(...
            remaining_isolated_objects));
        hypha_id = hypha_id_queaue(1);
    end
    hypha = (lbl_sk_st_isolated_hyphae == hypha_id);
    hypha_id_queaue(1) = [];
    connected = hypha | branch_mat;
    connected = bwareaopen(connected, 2);
    connected = connected & ~hypha;
    
    candidate_connectors = remaining_isolated_objects | connected;
    
    
    [xt, yt] = find(connected);
    % idx = sub2ind(img_dimension, yt, xt);
    disconnected = imfill(~candidate_connectors, [xt, yt], 8);
    
    candidate_connectors = candidate_connectors...
        & disconnected;
    candidate_connectors = candidate_connectors & ~connected;
    candidate_connectors = candidate_connectors & ~hypha;
    
    remaining_isolated_objects(hypha) = 0;
    
    connected_hypae_mat = lbl_sk_st_isolated_hyphae(candidate_connectors);
    
    connected_hyphae_id = unique(connected_hypae_mat);
    if ~isempty(connected_hyphae_id)
        hyphal_id_index = find(hypha_id==hyphal_labels);
        test = ismember(hyphal_labels, connected_hyphae_id);
        connected_hyphae_id_index = find(test);
        Adjacency_matrix(hyphal_id_index, connected_hyphae_id_index) = 1;
        % Symetric adjacency matrix
        Adjacency_matrix(connected_hyphae_id_index, hyphal_id_index) = 1;
    end
    hypha_id_queaue = [hypha_id_queaue; connected_hyphae_id];
end


G = graph(Adjacency_matrix,node_names);

end