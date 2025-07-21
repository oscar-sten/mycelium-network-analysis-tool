% A Ridge-Based Detection Algorithm with Filament Crossing
% Identification for 2D Mycelium Network Analysis


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
graph_representation = 1;
showLocalDensity = 0;
showLengthHistogram = 0;


theta = 0:15:360;
warning('off') % Suppress warnings


test_image = 's'; % Choose any letter in the italian alphabet.

suffix_img = '.jpg';


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
        sk_structure, [4, 4]);
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



%% Color based on ridge

if ridge_color_hyphae
    n_colors = 20;
    [~,E] = discretize(hyphal_median_ridge_values, n_colors);
    color_list = jet(n_colors);

    figure,
    imshow(I)
    hold on

    color_list = ["[0 0.4470 0.7410]", "[0.8500 0.3250 0.0980]",...
    "[0.9290 0.6940 0.1250]", "[0.4940 0.1840 0.5560]", ...
    "[0.4660 0.6740 0.1880]", "[0.3010 0.7450 0.9330]", ...
    "[0.6350 0.0780 0.1840]", "[0 0.6 0]", "[0 0 1]", ...
    "[1 0 0]", "[0 1 1]", "[1 1 1]"];

    for i=1:numel(hyphal_labels)
        hypha = (lbl_sk_st_isolated_hyphae==hyphal_labels(i));
        norm_median_ridge_value = hyphal_median_ridge_values(i);
        bin = find(E>norm_median_ridge_value);
        bin = bin(1)-1;
%         color = color_list(bin, :);

        color = color_list(mod(i, numel(color_list))+1);
        [X, Y] = find(hypha);
        scatter(Y, X, 'Marker', '.',...
            'MarkerFaceColor', color, ...
            'MarkerEdgeColor', color, 'HandleVisibility','off');

    end
    
%     colormap("jet")
%     colorbar("Ticks", [0, 1], "Ticklabels",...
%         {"Low", "High"}, "TickLabelInterpreter", "latex", "FontSize", 14)
end


%% Plot graph with centroid coordiniates


if graph_representation
    n_colors = 20;
    [~,E] = discretize(hyphal_median_ridge_values, n_colors);
    color_list = jet(n_colors);
    colors_markers = zeros(numel(hyphal_labels), 3);

    % Obtain centroid coordinates
    centroid_coords = zeros(numel(hyphal_labels), 2);

    figure,
    imshow(I)
    hold on
    for i=1:numel(hyphal_labels)
        label = hyphal_labels(i);
        hypha = (lbl_sk_st_isolated_hyphae==hyphal_labels(i));
        [hypha_x, hypha_y] = find(hypha);
        stats = regionprops(hypha);
        centroid = stats.Centroid;
        coords_tmp = [round(centroid(1)) round(centroid(2))];
        [K,D] = dsearchn([hypha_y, hypha_x], coords_tmp);
        centroid_coords(i, :) = [hypha_y(K), hypha_x(K)];
        norm_median_ridge_value = hyphal_median_ridge_values(i);
        bin = find(E>norm_median_ridge_value);
        bin = bin(1)-1;
        color = color_list(bin, :);
        colors_markers(i, :) = color;


    end
    % plot(centroid_coords(:, 1), centroid_coords(:, 2),'*')

    node_deg = degree(G);

    plot(G, 'XData', centroid_coords(:, 1), 'YData',...
        centroid_coords(:, 2),...
        'NodeColor', colors_markers, 'MarkerSize',...
        node_deg+1, 'EdgeColor',...
        'w', 'LineWidth', 1)
    colormap("jet")
    colorbar("Ticks", [0, 1], "Ticklabels",...
        {"Low", "High"}, "TickLabelInterpreter", "latex", "FontSize", 14)
    % title('Graph representation')
end


%% Local mycelium density
if showLocalDensity
    [local_density_map, global_density] = local_myclium_density(...
        sk_structure, mm_per_pixel);

    figure()
    imagesc(local_density_map);
    title('Density distribution')
    a = colorbar;
    % In Cardini et al. (2020) the mycelium density is defined
    % with the unit mm/mm^2.
    xlabel('[mm]')
    ylabel('[mm]')
    xticks([500, 1000, 1500, 2000, 2500])
    xticklabels(string(round([500, 1000, 1500, 2000, 2500]*mm_per_pixel, 2)))
    yticks([1920*(1/5), 1920*(2/5), 1920*(3/5), 1920*(4/5), 1920])
    yticklabels(string(round([1920*(1/5), 1920*(2/5), 1920*(3/5),...
        1920*(4/5), 1920]*mm_per_pixel, 2)))
    a.Label.String = 'mm/mm^2';
end

%% Length histogram
if showLengthHistogram
    for i=1:numel(hyphal_labels)
        hypha = (lbl_sk_st_isolated_hyphae==hyphal_labels(i));
        hy_lenght = compute_length(hypha, mm_per_pixel);
        hyphal_lengths(i) = hy_lenght;
    end
    figure()
    h = histogram(hyphal_lengths, 20, 'FaceColor', '[0 0.6 0]');
    ymax = max(h.BinCounts)*1.1;
    xlabel('Length [mm]')
    ylabel('Number of filaments')
    ylim([0 ymax])

end





