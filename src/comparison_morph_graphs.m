% Script for generating and comparing morphological graphs with and without
% overlaps. We consider the images s (#1 in the paper), o (#2 in the
% paper), b (#3 in the paper), p (#4 in the paper), r (#5 in the paper),
% and v (#6 in the paper).

clear all
addpath(genpath('fcn'))
addpath(genpath('..\data'))


%% Parameters

filter_sigma = 4;
D_T = 12;
patch_range = [10 10];
minimum_branch_length = 20;


theta_swich = 30; % Angle threshold of switching between functions for
% determining which is the secondary filament.

theta_T = 45; % Angle threshold for considering a pair of BPs likely
% to be a crossing.

theta_one = 90; % Angle threshold for detemining if a pair of BPs can be a
% crossing.


ridge_color_hyphae = 0;
morph_graph_representation = 0;
morph_graph_representation_w_overlaps = 0;
edge_weights_in_graph = 0;
plot_histograms = 1;

list_morph_edge_lengths = [];
list_morph_w_overlaps_edge_lengths = [];
list_allFilamentLengths = [];
list_node_ratios = [];
list_edge_ratios = [];

list_filament_degrees = [];

theta = 0:15:360;
warning('off') % Suppress warnings

suffix_img = '.jpg';

img_list = {'s', 'o', 'b', 'p', 'r', 'v', 'u', 'g'};

% img_list = {'r'};

for imgIdx = 1:numel(img_list)
    test_image = img_list{imgIdx};
    I = imread(strcat(test_image, suffix_img));

    %% Pre-analysis

    Igray = rgb2gray(I);

    if strcmp(test_image, 'h')
        mm_per_pixel = 0.0011;
    else
        [Igray, pixels_scale_bar, mm_scale_bar] = remove_scale_bar(...
            Igray);

        mm_per_pixel = mm_scale_bar/pixels_scale_bar;
    end

    nominal_hypal_widths_in_pixels = [floor(7.4/(1000*mm_per_pixel)) ...
        round(11.2/(1000*mm_per_pixel)) ceil(19/(1000*mm_per_pixel))];

    [polarity, vote] = sample_polarity(Igray);


    %% Mycelium detection

    [sk_structure, ridgeFilt] = mycelium_detection(Igray,...
        polarity, filter_sigma, theta, minimum_branch_length);

    sk_structure = bwskel(sk_structure);

    %% Remove disconnected components
    [sk_structure] = remove_disconnected_components(sk_structure);


    %% Analyse all BPs as branch points and obtaing direction vectors
    [branch_points, perfect_crossings_coords,  direction_vector_sets_3,...
        direction_vector_sets_4, connection_points] = get_BP_vectors(...
        sk_structure, ridgeFilt,...
        patch_range, theta_swich);


    %%  Obtain under passage coords for perfect crossings (Rare, Ad-hoc case)

    n_perfect_crossings = size(perfect_crossings_coords, 1);
    perfect_crossing_under_passage_coords = zeros(2*n_perfect_crossings, 2);

    for i=1:n_perfect_crossings
        bp = perfect_crossings_coords(i, :);
        under_passage_points = match_4_connected_BP(bp,...
            sk_structure, patch_range);
        perfect_crossing_under_passage_coords(2*i-1:2*i, :) =...
            under_passage_points;
    end


    %% Find tips
    tip_points = find_tip_points(sk_structure);



    %% Find sk neighbours

    [sk_neighbours, A] = find_sk_neigbours(sk_structure, branch_points,...
        tip_points);

    bp_keep = branch_points;
    tip_keep = tip_points;

    G = graph(A);


    %% Find possible crossing pairs
    [candidate_pairs, very_close_pairs, removed_pairs] = ...
        find_possible_crossing_pairs(...
        G, direction_vector_sets_3, branch_points, theta_one);


    %% Match very close pairs
    n_very_close_pairs = size(very_close_pairs, 1)/2;
    very_close_crossing_under_passage_coords = zeros(2*n_very_close_pairs, 2);

    for i=1:n_very_close_pairs
        bp = branch_points(very_close_pairs(2*i-1:2*i), :);
        under_passage_points = match_4_connected_BP(bp,...
            sk_structure, [3, 3]);
        very_close_crossing_under_passage_coords(2*i-1:2*i, :) =...
            under_passage_points;
    end

    very_close_pairs = unique(sort(very_close_pairs, 2), 'rows');

    %% Find likely crossing pairs
    [crossing_pairs, removed_pairs_step_2] = ...
        find_likely_crossing_pairs(G, direction_vector_sets_3,...
        candidate_pairs, branch_points, D_T, theta_T);



    %% Get coordinates

    [branch_points, connection_points, under_passage_points, crossing_centers] =...
        get_branch_and_crossing_coords(branch_points, connection_points,...
        perfect_crossings_coords,...
        crossing_pairs, perfect_crossing_under_passage_coords,...
        very_close_pairs, very_close_crossing_under_passage_coords);


    %% Separate filaments

    all_branches = connection_points;
    under_cross_coords = under_passage_points;
    [lbl_sk_st_isolated_hyphae, hyphal_labels, faulty_crossings] = ...
        isolate_filaments(sk_structure,...
        all_branches, under_cross_coords, 1);


    %% Compute median ridge of isolated hyphae
    hyphal_median_ridge_values = zeros(1, numel(hyphal_labels));

    ridgeImg = mat2gray(ridgeFilt);

    hyphal_lengths = zeros(numel(hyphal_labels), 1);

    for i=1:numel(hyphal_labels)
        idx = find(lbl_sk_st_isolated_hyphae==hyphal_labels(i));
        hy_lenght = numel(idx) * mm_per_pixel;
        hyphal_lengths(i) = hy_lenght;
        hy_ridge_resp = ridgeImg(idx);
        hyphal_median_ridge_values(i) = median(hy_ridge_resp);
    end


    %% Create graph object
    G = gen_graph(sk_structure, lbl_sk_st_isolated_hyphae, ...
        all_branches, under_cross_coords, hyphal_labels, ...
        hyphal_median_ridge_values);

    filament_degrees = degree(G);

    list_filament_degrees = [list_filament_degrees; filament_degrees];


    %% Original morphological graph without overlaps

    G_morph = graph(A);
    x_coords = [bp_keep(:, 1);tip_keep(:, 1)];
    y_coords = [bp_keep(:, 2); tip_keep(:, 2)];

    morph_node_coords = [x_coords, y_coords];


    G_morph.Nodes.Name = string((1:numel(x_coords))');
    G_morph.Nodes.Coordinates = [x_coords, y_coords];
    
    G_morph.Edges.Weight = G_morph.Edges.Weight*mm_per_pixel;

    %% Generate edge table for morphological graph with overlaps
    allEndNodes = {};
    allMiddleBPs = {};

    branchPointMap = zeros(size(lbl_sk_st_isolated_hyphae));
    bpMapIdx = sub2ind(size(lbl_sk_st_isolated_hyphae), branch_points(:, 2), ...
        branch_points(:, 1));
    branchPointMap(bpMapIdx) = 1;

    passagePointMap = zeros(size(lbl_sk_st_isolated_hyphae));
    passMapIdx = sub2ind(size(lbl_sk_st_isolated_hyphae),...
        under_cross_coords(:, 2), ...
        under_cross_coords(:, 1));
    passagePointMap(passMapIdx) = 1;

    allEdgesList = [];

    allFilamentLengths = [];

    for i=1:numel(hyphal_labels)
        hypha = (lbl_sk_st_isolated_hyphae==hyphal_labels(i));
        if sum(hypha(:)) > 3
            filament_endnode_coords = [];
            filament_endpoints = bwmorph(hypha, 'endpoints');
            [x, y] = find(filament_endpoints);
            filament_endpoint_coords = [y, x];
            % Check if tips (max 2)
            if any(ismember(filament_endpoint_coords, tip_points, 'rows'))
                idx = ismember(filament_endpoint_coords, tip_points, 'rows');
                if sum(idx) > 2
                    % In some cases, multiple endpoints are detected in at the end
                    % of the filament, out of these all but one has a degree of 0
                    % in the morphological graph, and can thus be identified and
                    % removed
                    coords_tmp = filament_endpoint_coords(idx, :);
                    idx_tmp = ismember(morph_node_coords, coords_tmp, 'rows');
                    cand_tip_idx = find(idx_tmp);
                    deg_morph = degree(G_morph, cand_tip_idx);
                    test = cand_tip_idx(deg_morph==0);
                    test_coords = morph_node_coords(test, :);
                    idx(ismember(filament_endpoint_coords,...
                        test_coords, "rows")) = 0;
                end
                filament_endnode_coords = [filament_endnode_coords;...
                    filament_endpoint_coords(idx, :)];
            end

            % Check if any connection points (max 2) (max distance sqrt(2)
            if size(filament_endnode_coords, 1) < 2
                [k, dist] = dsearchn(connection_points, [y, x]);
                close_connection_pts_idx = (dist <= sqrt(2));
                close_connection_pts = k(close_connection_pts_idx);
                filament_endnode_coords = [filament_endnode_coords;...
                    branch_points(close_connection_pts, :)];
                filament_endnode_connection_pts = connection_points(...
                    close_connection_pts, :);
            else
                filament_endnode_connection_pts = [];
            end

            allEndNodes{end+1} = filament_endnode_coords;

            if size(filament_endnode_coords, 1) < 2
                error("Too few end points")
            elseif size(filament_endnode_coords, 1) > 2
                error("Too many end points")
            end

            % finde the branch points on the filament
            bpsOnFilament = hypha & branchPointMap;
            [x, y] = find(bpsOnFilament);
            bpsOnFilament_coords = [y, x];
            allMiddleBPs{end+1} = bpsOnFilament_coords;

            % Add endnodes to filament
            if ~isempty(filament_endnode_connection_pts)
                filament_endnode_idx = sub2ind(size(hypha),...
                    [filament_endnode_coords(:, 2);...
                    filament_endnode_connection_pts(:, 2)], ...
                    [filament_endnode_coords(:, 1);...
                    filament_endnode_connection_pts(:, 1)]);
                hypha(filament_endnode_idx) = 1;
            end

            % It is important that no inappropriate shortcuts are taken (for
            % instance, if a filament goes in a loop, the an end point can
            % be closer to the other end point than to a BP in the middle,
            % therefore, for each "under" passage in the filament,
            % the connectivity is restored by drawing a straight line between
            % the passage points. Then the geodesic distance transform is used
            % to deterimine the edge weights and the topological order.

            % First we isolate the passage points for the particular filament.
            currentHyphaWPassagePts = (hypha | passagePointMap);
            currentHyphaWPassagePts = bwareaopen(currentHyphaWPassagePts, 3);
            currentPassagePts = (currentHyphaWPassagePts & passagePointMap);
            [x, y] = find(currentPassagePts);
            currentPassageCoords = [y, x];


            % Maybe, assure that the pairs of points are correct?
            % Also add criterion to only build the bridge if it's needed.
            % Fill the gaps between the under passage points in the filament.
            for ii=1:2:size(currentPassageCoords, 1)-1
                points = currentPassageCoords(ii:ii+1, :);
                % Connect two points
                nPoints = max(abs(points(1, 1) - points(2, 1)), ...
                    abs(points(1, 2) - points(2, 2)))+1;

                rIndex = round(linspace(points(1, 2), ...
                    points(2, 2), nPoints));  % Row indices

                cIndex = round(linspace(points(1, 1),...
                    points(2, 1), nPoints));  % Column indices

                index = sub2ind(size(hypha), rIndex, cIndex);
                hypha(index) = 1;
            end

            % Compute geodesic distance transform with first end point as seed.
            % Quasi-eucliden counts diagonal steps as sqrt(2).
            geodesicDist = bwdistgeodesic(hypha, filament_endnode_coords(1, 1),...
                filament_endnode_coords(1, 2), 'quasi-euclidean');

            current_filament_length = geodesicDist(...
                filament_endnode_coords(2, 2), filament_endnode_coords(2, 1));

            allFilamentLengths = [allFilamentLengths; current_filament_length];

            % Build local subgraph staring with seed as starting point
            seed_node = [filament_endnode_coords(1, 1),...
                filament_endnode_coords(1, 2)];
            other_nodes = [bpsOnFilament_coords;...
                filament_endnode_coords(2, 1),...
                filament_endnode_coords(2, 2)];
            other_nodes_idx = sub2ind(size(hypha), other_nodes(:, 2),...
                other_nodes(:, 1));
            other_nodes_dists = geodesicDist(other_nodes_idx);
            edge_list = zeros(size(other_nodes_dists, 1), 5);
            starting_node = seed_node;
            for ii=1:size(other_nodes_dists, 1)
                [dist, idx] = min(other_nodes_dists);
                end_node = other_nodes(idx, :);
                edge_list(ii, 1:2) = starting_node;
                edge_list(ii, 3:4) = end_node;
                edge_list(ii, 5) = dist;
                starting_node = end_node;
                other_nodes_dists = other_nodes_dists-dist;
                other_nodes_dists(other_nodes_dists<=0) = Inf;
            end

            allEdgesList = [allEdgesList; edge_list];
        end

    end

    %% Build new graph

    edge_list = allEdgesList;

    % Generate node table
    allNodeCoords = unique([edge_list(:, 1:2); edge_list(:, 3:4)], "rows");

    nodeIDs = string((1:size(allNodeCoords, 1))');

    NodeTable = table(nodeIDs, allNodeCoords,'VariableNames',...
        {'Name' 'Coordinates'});

    % Generate edge table
    source_list = zeros(size(edge_list, 1), 1);
    target_list = zeros(size(edge_list, 1), 1);
    weight_list = zeros(size(edge_list, 1), 1);

    for i=1:size(edge_list, 1)
        % Identify and store the node indices of source and target node in each
        % edge.
        current_edge = edge_list(i, :);
        current_source = find(ismember(allNodeCoords, current_edge(1:2),...
            'rows'));
        current_target = find(ismember(allNodeCoords, current_edge(3:4),...
            'rows'));
        current_weight = current_edge(5);
        % Store
        source_list(i) = current_source;
        target_list(i) = current_target;
        weight_list(i) = current_weight;
    end

    weight_list = weight_list*mm_per_pixel;

    EdgeTable = table([source_list target_list], weight_list, ...
        'VariableNames',{'EndNodes' 'Weight'});


    G_morph_w_overlaps = graph(EdgeTable, NodeTable);


    %% Apply corrections
    allEdgesWeights = G_morph_w_overlaps.Edges.Weight;
    invalidEdgesIdx = find((allEdgesWeights==0) | (allEdgesWeights==Inf));
    G_morph_w_overlaps = rmedge(G_morph_w_overlaps, invalidEdgesIdx);
    G_morph_w_overlaps = removeBridgeNodes(G_morph_w_overlaps);

    % Test if node degrees are correct
    if sum(degree(G_morph_w_overlaps, G_morph_w_overlaps.Nodes.Name)==2)~=0 ||...
            sum(degree(G_morph_w_overlaps, G_morph_w_overlaps.Nodes.Name)>3)
        error("Wrong degree")
    end


    %% Plotting
    if ridge_color_hyphae
        n_colors = 20;
        [~,E] = discretize(hyphal_median_ridge_values, n_colors);
        color_list = jet(n_colors);

        figure,
        imshow(I)
        hold on

        for i=1:numel(hyphal_labels)
            hypha = (lbl_sk_st_isolated_hyphae==hyphal_labels(i));
            norm_median_ridge_value = hyphal_median_ridge_values(i);
            bin = find(E>norm_median_ridge_value);
            bin = bin(1)-1;
            color = color_list(bin, :);
            [X, Y] = find(hypha);
            scatter(Y, X, 'Marker', '.',...
                'MarkerFaceColor', color, ...
                'MarkerEdgeColor', color, 'HandleVisibility','off');

        end

        plot(under_cross_coords(:, 1), under_cross_coords(:, 2), 'go')

    end

    if morph_graph_representation
        x_coords = G_morph.Nodes.Coordinates(:, 1);
        y_coords = G_morph.Nodes.Coordinates(:, 2);
        figure,
        imshow(I)
        hold on
        if ~edge_weights_in_graph
            plot(G_morph, 'XData', x_coords, 'YData', y_coords, 'LineStyle',...
                '-','NodeLabel',{},...
                'EdgeColor', 'b', 'NodeColor', 'r', 'LineWidth', 1.5, ...
                'MarkerSize', 2)
        else
            plot(G_morph, 'XData', x_coords, 'YData', y_coords, 'LineStyle',...
                '-','NodeLabel',{},...
                'EdgeColor', 'b', 'NodeColor', 'r', 'LineWidth', 1.5, ...
                'MarkerSize', 2, 'EdgeLabel', round(G_morph.Edges.Weight, 2))
        end
    end


    if morph_graph_representation_w_overlaps
        x_coords = G_morph_w_overlaps.Nodes.Coordinates(:, 1);
        y_coords = G_morph_w_overlaps.Nodes.Coordinates(:, 2);
        figure,
        imshow(I)
        hold on
        if ~edge_weights_in_graph
            plot(G_morph_w_overlaps, 'XData', x_coords, 'YData', y_coords,...
                'LineStyle', '-','NodeLabel',{},...
                'EdgeColor', 'b', 'NodeColor', 'r', 'LineWidth', 1.5, ...
                'MarkerSize', 2)
        else
            plot(G_morph_w_overlaps, 'XData', x_coords, 'YData', y_coords,...
                'LineStyle', '-','NodeLabel',{},...
                'EdgeColor', 'b', 'NodeColor', 'r', 'LineWidth', 1.5, ...
                'MarkerSize', 2, 'EdgeLabel',...
                round(G_morph_w_overlaps.Edges.Weight, 2))
        end
    end
    
    %% Comparisons betwee repressentations
    % Node numbers
    nodeDegMorph = degree(G_morph);
    nodeDegMorphWOverlaps = degree(G_morph_w_overlaps);
    edgeLengthsMorph = G_morph.Edges.Weight;
    edgeLengthsMorphWOverlaps = G_morph_w_overlaps.Edges.Weight;
    allFilamentLengths = allFilamentLengths*mm_per_pixel;

    % Stack results
    list_morph_edge_lengths = [list_morph_edge_lengths; edgeLengthsMorph];
    list_morph_w_overlaps_edge_lengths = ...
        [list_morph_w_overlaps_edge_lengths; edgeLengthsMorphWOverlaps];
    list_allFilamentLengths = [list_allFilamentLengths; allFilamentLengths];

    %% Print stats

    disp(strcat("------------- Image: ", test_image, " ---------------"))
    N_g = sum(nodeDegMorph==3)+sum(nodeDegMorph==2);
    N = sum(nodeDegMorphWOverlaps==3);

    disp(strcat("Number of branch nodes in original morphological graph: ",...
        string(N_g)))

    disp(strcat("Number of branch nodes in improved graph: ",...
        string(N)))

    ratio_nodes = 1-(N/N_g);

    disp(strcat(...
        "Ratio biological nodes over geometrical (Dikec et al. 2020): ", ...
        string(round(ratio_nodes,2))))

    disp(strcat("Number of edges in original morphological graph: ",...
        string(numel(edgeLengthsMorph))))
    disp(strcat("Number of edges in morphological graph with overlaps: ",...
        string(numel(edgeLengthsMorphWOverlaps))))

    disp(strcat("Ratio number of edges: ",...
        string(round(numel(edgeLengthsMorphWOverlaps)/...
        numel(edgeLengthsMorph), 2))))


    ratio_edges = 1-(numel(edgeLengthsMorphWOverlaps)/...
        numel(edgeLengthsMorph));



    list_node_ratios = [list_node_ratios; ratio_nodes];
    list_edge_ratios = [list_edge_ratios; ratio_edges];


end

%% Plot histograms
if plot_histograms
    figure,
    histogram(list_morph_edge_lengths, 'BinWidth', 0.05, "FaceColor",...
        '[0, 0.7, 0]')
    %title('Edge lengths in original morphological graph')
    xlabel("Lenght [mm]")
    ylabel("Number of edges")
    ylim([0, 750])

    figure,
    histogram(list_morph_w_overlaps_edge_lengths, 'BinWidth', 0.05,...
        "FaceColor", '[0, 0.7, 0]')
    %title('Edge lengths in morphological graph with overlaps')
    xlabel('Lenght [mm]')
    ylabel("Number of edges")
    ylim([0, 750])

%     figure,
%     histogram(list_allFilamentLengths)
%     title('Lengths of segmented filaments')
%     xlabel('Lenght [mm]')
%% djsk

    figure,
    % k_min = 1
    hist = histogram(list_filament_degrees, "FaceColor", '[0, 0.7, 0]');
%     hold on
%     Y = hist.Values';
%     Y = Y(2:21);
% 
%     X = (1:20)';
%     [f, gof] = fit(X, Y, 'power1');
%     x_range = linspace(1, 20);
%     [alpha, xmin, L]=plfit(list_filament_degrees(list_filament_degrees>0), ...
%         'xmin', 1);
% %     [p,gof]=plpva(list_filament_degrees(list_filament_degrees>0), 5);
% %     % xmin=5 highest p=value.
% 
% 
%     plot(x_range, 307*(x_range.^-alpha), 'r', 'LineWidth', 2)
    ylabel("Number of filaments")
    xlabel("Degree")
end


% %% dksl
% 
% divisor(1) = 1^f.b + 2^f.b;
% for i=2:20
%     divisor(i) = divisor(i-1) + (i+1)^f.b;
% end
% 
% pr = (1./divisor)'.*(X.^f.b);
% 
% Y_pdf = Y/sum(Y);
% 
% [alpha, xmin, L]=plfit(list_filament_degrees(list_filament_degrees>0));
% 

%% Print final stats

disp("------------------ Summary ------------------")
disp(strcat(...
    "Mean of Dikec et al. (2020) r-value: ",...
    string(round(mean(list_node_ratios), 2))))
disp(strcat(...
    "Standard deviation of Dikec et al. (2020) r-value: ",...
    string(round(std(list_node_ratios), 2))))

disp(strcat("Mean ratio edge number: ",...
    string(round(mean(list_edge_ratios), 2))))





