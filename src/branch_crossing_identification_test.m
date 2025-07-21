% Branch and crossing identification test

clear all
addpath(genpath('fcn'))
addpath(genpath('..\data\lbl_intersections_test\'))


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



%% Test


suffix_annotation = '_intersections_label.mat';


GT = load(strcat(test_image, suffix_annotation));
intersection_labels = GT.gTruth.LabelData;

figure,
imshow(I)
hold on
spy(sk_structure, 'b')
xlabel('')

% Quantify accuracy

% Rectangle coordinates are in fromat [x y w h] where x and y define the
% upper left corner and w is width and h is height.

% Branches
all_branches = branch_points;
identified_branches = all_branches;
correctly_identified_branches = [];
false_negative_branches = [];
branchBoxes = table2array(intersection_labels(:, 'Branch'));
branchBoxes = branchBoxes{1};

n_branches = size(branchBoxes, 1);

for i=1:n_branches
    currentBox = branchBoxes(i, :);
    x_interval = [currentBox(1), currentBox(1)+currentBox(3)];
    y_interval = [currentBox(2), currentBox(2)+currentBox(4)];
    idx = find((identified_branches(:, 1)>x_interval(1)) & ...
        identified_branches(:, 1)<x_interval(2) & ...
        identified_branches(:, 2)>y_interval(1) & ...
        identified_branches(:, 2)<y_interval(2));
    if ~isempty(idx)
        correct_branch = identified_branches(idx, :);
        correctly_identified_branches = [correctly_identified_branches;...
            correct_branch];
        identified_branches(idx, :) = [];
    else
        false_negative_branches = [false_negative_branches;...
            [mean(x_interval), mean(y_interval)]];
    end
    rectangle('Position', branchBoxes(i, :), 'EdgeColor','[0, 0.6, 0]');
end

plot(correctly_identified_branches(:, 1),...
    correctly_identified_branches(:, 2), 'o', 'Color', 'k')

if ~isempty(false_negative_branches)
    plot(false_negative_branches(:, 1), false_negative_branches(:, 2), '^', ...
        'Color', '[0.9290 0.6940 0.1250]', 'LineWidth', 2)
end

if ~isempty(identified_branches)
    plot(identified_branches(:, 1), identified_branches(:, 2), '^', ...
        'Color', 'red', 'LineWidth', 2)
end


% Crossings
identified_crossings = crossing_centers;
correctly_identified_crossings = [];
false_negative_crossings = [];
crossingBoxes = table2array(intersection_labels(:, "Crossing"));
crossingBoxes = crossingBoxes{1};
n_crossings = size(crossingBoxes, 1);

for i=1:n_crossings
    currentBox = crossingBoxes(i, :);
    x_interval = [currentBox(1), currentBox(1)+currentBox(3)];
    y_interval = [currentBox(2), currentBox(2)+currentBox(4)];
    idx = find((identified_crossings(:, 1)>x_interval(1)) & ...
        identified_crossings(:, 1)<x_interval(2) & ...
        identified_crossings(:, 2)>y_interval(1) & ...
        identified_crossings(:, 2)<y_interval(2));
    if ~isempty(idx)
        correct_crossing = identified_crossings(idx, :);
        correctly_identified_crossings = [correctly_identified_crossings;...
            correct_crossing];
        identified_crossings(idx, :) = [];
    else
        false_negative_crossings = [false_negative_crossings;...
            [mean(x_interval), mean(y_interval)]];
    end
    rectangle('Position', crossingBoxes(i, :), 'EdgeColor','b');
end


plot(correctly_identified_crossings(:, 1),...
    correctly_identified_crossings(:, 2), 'o', 'Color', 'k')

if ~isempty(false_negative_crossings)
    plot(false_negative_crossings(:, 1), false_negative_crossings(:, 2),...
        's', 'Color', '[0.9290 0.6940 0.1250]', 'LineWidth', 2)
end

if ~isempty(identified_crossings)
    plot(identified_crossings(:, 1), identified_crossings(:, 2), ...
        's', 'Color', 'red', 'LineWidth', 2)
end


n_tp_branch = size(correctly_identified_branches, 1);
n_fp_branch = size(identified_branches, 1);
n_fn_branch = size(false_negative_branches, 1);

n_tp_crossing = size(correctly_identified_crossings, 1);
n_fp_crossing = size(identified_crossings, 1);
n_fn_crossing = size(false_negative_crossings, 1);

figure()
groups = categorical(["Branches", "Crossings"]);
b = bar(groups, [n_tp_branch, n_fp_branch, n_fn_branch;...
    n_tp_crossing, n_fp_crossing, n_fn_crossing], 'grouped');
%% jfskl
b(1).FaceColor = [0, 0.6, 0];
b(2).FaceColor = 'red';
b(3).FaceColor = [0.9290 0.6940 0.1250];

% Put text above bars
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips3 = b(3).XEndPoints;
ytips3 = b(3).YEndPoints;
labels3 = string(b(3).YData);
text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')


legend('True positive', 'False positive', 'False negative')
%% kdsl
% Precision
precision_branches = n_tp_branch/(n_tp_branch + n_fp_branch);
precision_crossings = n_tp_crossing/(n_tp_crossing + n_fp_crossing);

% Recall
recall_branches = n_tp_branch/(n_tp_branch + n_fn_branch);
recall_crossings = n_tp_crossing/(n_tp_crossing + n_fn_crossing);

% F1 score
f1_branches = (2*n_tp_branch)/(2*n_tp_branch + n_fn_branch + n_fp_branch);
f1_crossings = (2*n_tp_crossing)/(2*n_tp_crossing +...
    n_fn_crossing + n_fp_crossing);

% Print stats

disp(strcat("Precision in identifying branches: ", ...
    string(round(precision_branches*100, 2)), "%"))

disp(strcat("Recall in identifying branches: ", ...
    string(round(recall_branches*100, 2)), "%"))

disp(strcat("F1 score on identifying branches: ", ...
    string(round(f1_branches, 2))))

disp(strcat("Precision in identifying crossings: ", ...
    string(round(precision_crossings*100, 2)), "%"))

disp(strcat("Recall in identifying crossings: ", ...
    string(round(recall_crossings*100, 2)), "%"))

disp(strcat("F1 score on identifying crossings: ", ...
    string(round(f1_crossings, 2))))



