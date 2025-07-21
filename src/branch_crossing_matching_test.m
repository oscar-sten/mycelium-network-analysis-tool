% Matching test of branches and crossings

clear all
addpath(genpath('fcn'))
addpath(genpath('..\data\Cardini_et_al_images\'))
addpath(genpath('..\data\lbl_matching_test\'))


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

theta = 0:15:360;
warning('off') % Suppress warnings

% OBS: ANNOTATION ERROR
test_image = 's';

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



%% Test matching accuracy 
patch_range_single = 15;
patch_range_pair = 20;
patch_range_multi = 30;
GT = load('s_lbl_test.mat');
im_gt = imread(GT.gTruth.LabelData.PixelLabelData{1});
% Color based on ridge

n_colors = 20;
[~,E] = discretize(hyphal_median_ridge_values, n_colors);
color_list = jet(n_colors);

figure,
imshow(I)
hold on

% Ground truth labels
spy(im_gt==1, 'r') % Main hypha

spy(im_gt==2, 'm') % Branches

spy(im_gt==3, 'g') % Under passages

for i=1:numel(hyphal_labels)
    hypha = (lbl_sk_st_isolated_hyphae==hyphal_labels(i));
    norm_median_ridge_value = hyphal_median_ridge_values(i);
    %     color = [0, 1-norm_median_ridge_value, 0];
    bin = find(E>norm_median_ridge_value);
    bin = bin(1)-1;
    color = color_list(bin, :);
    [X, Y] = find(hypha);
    scatter(Y, X, 'Marker', '.',...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', color);

end


plot(all_branches(:, 1), all_branches(:, 2), 'mo')

plot(under_cross_coords(:, 1), under_cross_coords(:, 2), 'go')


xlabel('')

%%  Quantify error of branches

[n_branches, ~] = size(all_branches);
[n,m] = size(Igray);

correct_branches = 0;
incorrect_branches = 0;
correct_bp_set = [];
incorrect_bp_set = [];

for i = 1:n_branches
    bp = all_branches(i, :);
    if im_gt(bp(2), bp(1))~=0
        % Define patch indices
        x_range = bp(1)-patch_range_single:bp(1)+patch_range_single;
        y_range = bp(2)-patch_range_single:bp(2)+patch_range_single;

        % Remove invalid indices
        x_range(x_range<1)=[];
        x_range(x_range>m)=[];
        y_range(y_range<1)=[];
        y_range(y_range>n)=[];

        % Extract patches
        gt_patch = im_gt(y_range, x_range);
        lbl_patch = lbl_sk_st_isolated_hyphae(y_range, x_range);

        % Remove hyphae not connected to the branch point from lbl patch
        in_patch = bp - [x_range(1) y_range(1)] + [1 1];

        in_patch_bw = zeros(size(lbl_patch));
        xt = in_patch(:, 1);
        yt = in_patch(:, 2);
        idx = sub2ind(size(lbl_patch), yt, xt);
        in_patch_bw(idx) = 1;

        % Set bp as 1 to connect it
        lbl_patch(idx) = 1;

        disconnected = imfill(~lbl_patch, [in_patch(:, 2) in_patch(:, 1)], 8); % Test if correct

        % Set bp as zero to remove connection again
        lbl_patch(idx) = 0;

        lbl_patch(~disconnected) = 0;

        % GT patch parts
        main_part = (gt_patch==1);
        brach_part = (gt_patch==2);

        % Identify and count overlapping pixels between annotation and
        % different labels
        labels = unique(lbl_patch);
        main_pxl = zeros(1, numel(labels)-1);
        branch_pxl = zeros(1, numel(labels)-1);
        for k=2:numel(labels)
            hypha = lbl_patch==labels(k);
            main_corresp_patch = (main_part & hypha);
            brach_corresp_patch = (brach_part & hypha);
            main_pxl(k-1) = sum(main_corresp_patch(:));
            branch_pxl(k-1) = sum(brach_corresp_patch(:));
        end

        % Test if the pixel overlappings are significantly diffrent from
        % each other.
        test_main = abs(main_pxl(1)-main_pxl(2))/...
            (max(main_pxl(1), main_pxl(2)));
        test_branch = abs(branch_pxl(1)-branch_pxl(2))/...
            (max(branch_pxl(1), branch_pxl(2)));
        if (test_main > 0.5) && (test_branch > 0.5)
            correct_branches = correct_branches + 1;
            correct_bp_set = [correct_bp_set; bp];
        else
            incorrect_branches = incorrect_branches + 1;
            incorrect_bp_set = [incorrect_bp_set; bp];
        end

    end
end

if ~isempty(correct_bp_set)
    plot(correct_bp_set(:, 1), correct_bp_set(:, 2), '^', 'Color',...
        '[0, 0.6, 0]', 'LineWidth', 2)
end

if ~isempty(incorrect_bp_set)
    plot(incorrect_bp_set(:, 1), incorrect_bp_set(:, 2), '^', 'Color',...
        'r', 'LineWidth', 2)
    t=1;
end


% Quantify error of crossings

n_crossings = size(crossing_centers, 1);

correct_crossings = 0;
incorrect_crossings = 0;

correct_crossing_set = [];
incorrect_crossing_set = [];

for i = 1:n_crossings

    crossing_center = ceil(crossing_centers(i, :));
    if im_gt(crossing_center(2), crossing_center(1))~=0

        % Define patch indices
        x_range = crossing_center(1)-...
            patch_range_single:crossing_center(1)+patch_range_single;
        y_range = crossing_center(2)-...
            patch_range_single:crossing_center(2)+patch_range_single;

        % Remove invalid indices
        x_range(x_range<1)=[];
        x_range(x_range>m)=[];
        y_range(y_range<1)=[];
        y_range(y_range>n)=[];

        % Extract patches
        gt_patch = im_gt(y_range, x_range);
        lbl_patch = lbl_sk_st_isolated_hyphae(y_range, x_range);

        % GT patch parts
        main_part = (gt_patch==1);
        under_passing_part = (gt_patch==3);

        % Identify and count overlapping pixels between annotation and
        % different labels
        labels = unique(lbl_patch);
        main_pxl = zeros(1, numel(labels)-1);
        under_pass_pxl = zeros(1, numel(labels)-1);
        for k=2:numel(labels)
            hypha = lbl_patch==labels(k);
            main_corresp_patch = (main_part & hypha);
            under_passing_patch = (under_passing_part & hypha);
            main_pxl(k-1) = sum(main_corresp_patch(:));
            under_pass_pxl(k-1) = sum(under_passing_patch(:));
        end
        % Test if the pixel overlappings are significantly diffrent from
        % each other.
        % Replace with something statistically rigorous.
        test_main = abs(main_pxl(1)-main_pxl(2))/...
            (max(main_pxl(1), main_pxl(2)));
        test_under_pass = abs(under_pass_pxl(1)-under_pass_pxl(2))/...
            (max(under_pass_pxl(1), under_pass_pxl(2)));
        if (test_main > 0.5) && (test_under_pass > 0.5)
            correct_crossings = correct_crossings + 1;
            correct_crossing_set = [correct_crossing_set;...
                crossing_center];
        else
            incorrect_crossings = incorrect_crossings + 1;
            incorrect_crossing_set = [incorrect_crossing_set;...
                crossing_center];
        end
    end
end

if ~isempty(correct_crossing_set)
    plot(correct_crossing_set(:, 1), ...
        correct_crossing_set(:, 2), 'square', 'Color',...
        '[0, 0.6, 0]', 'LineWidth', 2)
end

if ~isempty(incorrect_crossing_set)
    plot(incorrect_crossing_set(:, 1), ...
        incorrect_crossing_set(:, 2), 'square', 'Color',...
        'r', 'LineWidth', 2)
end


%% Bar chart
figure()
x = categorical({'Branches', 'Crossings'});
y = [correct_branches incorrect_branches;...
    correct_crossings incorrect_crossings];

b = bar(x, y, 'stacked', 'FaceColor','flat');
b(1).CData = [0 0.6 0];
b(2).CData = [0.8 0 0];


%% Print

% For image S: 160 (152) branches and 20 (17) crossings.


disp("------------------------ Results ---------------------------")
disp(strcat("Number of evaluated branches: ", ...
    string(correct_branches+incorrect_branches)))

disp(strcat("Number of correct branches: ", ...
    string(correct_branches)))

disp(strcat("Number of incorrect branches: ", ...
    string(incorrect_branches)))

disp(strcat("Accuracy of matching for branches: ", ...
    string(...
    round(correct_branches/(correct_branches+incorrect_branches)*100,...
    1)),"%"))

disp(strcat("Number of evaluated crossings: ", ...
    string(correct_crossings + incorrect_crossings)))

disp(strcat("Number of correct crossings: ", ...
    string(correct_crossings)))

disp(strcat("Number of incorrect crossings: ", ...
    string(incorrect_crossings)))

disp(strcat("Accuracy of matching for crossings: ",...
    string(...
    round(correct_crossings/(correct_crossings+incorrect_crossings)*100,...
    1)), "%"))














