% Test the segmentation performance on the images created by Cardini and
% annotated by Sten.

% Obs: RAM requirement for cost matrix ~ 57692x61104 (26.3GB)
% This script is very computationally demanding. 
% We recomend atleast 128 GB of RAM and Intel i9 processor or better.

clear all
%close all
% clc
addpath(genpath('fcn')) % Add functions to path

data_path_GT = '..\data\ground_truth_skeletons_Cardini_et_al_images\';
addpath(genpath(data_path_GT))

data_path_img = '..\data\Cardini_et_al_images\';
addpath(genpath(data_path_img))

result_save_path = '..\result_tables\';
addpath(genpath(result_save_path))


image_list = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'l',...
    'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v'];

image_prefix = '';
image_suffux = '.jpg';

GT_prefix = 'GT_skeleton_';
GT_suffix = '.png';


precision_list = [];

recall_list = [];

F1_list = [];

MMC_list = [];

analyzed_images = [];

n_TP_List = [];
n_TN_List = [];
n_FP_List = [];
n_FN_List = [];

%% Parameters
tolerance = 8.4853;
minimum_branch_length = 20;
filter_sigma = 4;
result_decimals = 3;
plot = 0;
save = 1;


disp(strcat("Tolerance: ", string(tolerance)))


for i = 1:numel(image_list)

    I = imread(strcat(image_prefix, image_list(i), image_suffux));
    
    disp(strcat("Processing image: ", image_list(i)));
    % Image name for table
    analyzed_images = [analyzed_images; strcat(image_prefix,...
        image_list(i))];

    GT = imread(strcat(GT_prefix, image_list(i), GT_suffix));

    % Dilate GT for calculating precision, dilate sk_structure for
    % caluculating recall.

    GT_bin = imbinarize(mat2gray(GT));

    %     GT_bin = ~GT_bin;

    Igray = rgb2gray(I);
    if strcmp(image_list(i), 'h')
        mm_per_pixel = 0.0011;
    else
        [Igray, pixels_scale_bar, mm_scale_bar] = remove_scale_bar(...
            Igray);

        mm_per_pixel = mm_scale_bar/pixels_scale_bar;
    end

    [polarity, vote] = sample_polarity(Igray);

    theta = [0:15:360];
    
    [sk_structure, ridgeFilt] = mycelium_detection(Igray,...
        polarity, filter_sigma, theta, minimum_branch_length);


    %% Remove disconnected components
    [sk_structure] = remove_disconnected_components(sk_structure);

    
 % Find the indices of non-zero elements
    [rowA, colA] = find(sk_structure);
    [rowB, colB] = find(GT_bin);
    
    disp("Building cost matrix")
    % Initialize the cost matrix
    costMatrix = zeros(length(rowA), length(rowB));

    % Calculate the Euclidean distance between all pairs of non-zero elements
    for k = 1:length(rowA)
        for j = 1:length(rowB)
            dist = sqrt((rowA(k)-rowB(j))^2 + (colA(k)-colB(j))^2);
            costMatrix(k,j) = dist;
        end
    end
    
    disp("Solving optimal matching problem")
    % Solve assignment problem
    [M,uR,uC] = matchpairs(costMatrix, tolerance/2);

    [n, m] = size(sk_structure);

    % Quantify
    TP = size(M, 1);
    FN = size(uC, 1);
    FP = size(uR, 1);
    TN = n*m - TP - FN - FP;

    % Compute stats
    precision = TP/(TP+FP);
    recall = TP/(TP+FN);
    F1 = (2*TP)/(2*TP+FP+FN);
    MMC = (TP*TN - FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)));

    % Store
    precision_list = [precision_list; precision];
    recall_list = [recall_list; recall];
    F1_list = [F1_list; F1];
    MMC_list = [MMC_list; MMC];

    % Store TP, FN, FP, and TN
    n_TP_List = [n_TP_List; TP];
    n_TN_List = [n_TN_List; TN];
    n_FP_List = [n_FP_List; FP];
    n_FN_List = [n_FN_List; FN];
    
    disp("Done")

end


%% Summarize results in a table

% Compute averages
average_precision = mean(precision_list);
average_recall = mean(recall_list);
average_F1 = mean(F1_list);
average_MMC = mean(MMC_list);
average_TP = mean(n_TP_List);
average_FN = mean(n_FN_List);
average_FP = mean(n_FP_List);
average_TN = mean(n_TN_List);

% Add averages to tables
analyzed_images = [analyzed_images; "Average"];
precision_list = [precision_list; average_precision];
recall_list = [recall_list; average_recall];
F1_list = [F1_list; average_F1];
MMC_list = [MMC_list; average_MMC];
n_TP_List = [n_TP_List; average_TP];
n_FN_List = [n_FN_List; average_FN];
n_FP_List = [n_FP_List; average_FP];
n_TN_List = [n_TN_List; average_TN];


% Compute medians 
median_precision = median(precision_list);
median_recall = median(recall_list);
median_F1 = median(F1_list);
median_MMC = median(MMC_list);
median_TP = median(n_TP_List);
median_FN = median(n_FN_List);
median_FP = median(n_FP_List);
median_TN = median(n_TN_List);


% Add medians to tables
analyzed_images = [analyzed_images; "Median"];
precision_list = [precision_list; median_precision];
recall_list = [recall_list; median_recall];
F1_list = [F1_list; median_F1];
MMC_list = [MMC_list; median_MMC];
n_TP_List = [n_TP_List; median_TP];
n_FN_List = [n_FN_List; median_FN];
n_FP_List = [n_FP_List; median_FP];
n_TN_List = [n_TN_List; median_TN];

% Add sums
analyzed_images = [analyzed_images; "Sum"];
sum_TP = sum(n_TP_List(1:numel(image_list)));
sum_FP = sum(n_FP_List(1:numel(image_list)));
sum_FN = sum(n_FN_List(1:numel(image_list)));
sum_TN = sum(n_TN_List(1:numel(image_list)));
n_TP_List = [n_TP_List; sum_TP];
n_FN_List = [n_FN_List; sum_FN];
n_FP_List = [n_FP_List; sum_FP];
n_TN_List = [n_TN_List; sum_TN];

sum_precision = sum_TP/(sum_TP + sum_FP);
sum_recall = sum_TP/(sum_TP + sum_FN);
sum_F1 = (2*sum_TP)/(2*sum_TP + sum_FN + sum_FP);
sum_MMC = (sum_TP*sum_TN -...
        sum_FP*sum_FN)/(sqrt((sum_TP+sum_FP)*...
        (sum_TP+sum_FN)*(sum_TN+sum_FP)*...
        (sum_TN+sum_FN)));

precision_list = [precision_list; sum_precision];
recall_list = [recall_list; sum_recall];
F1_list = [F1_list; sum_F1];
MMC_list = [MMC_list; sum_MMC];


% Round to decired number of decimals
precision_list = round(precision_list, result_decimals);
recall_list = round(recall_list, result_decimals);
F1_list = round(F1_list, result_decimals);
MMC_list = round(MMC_list, result_decimals);


% Create table
varNames = ["Image", "Precision", "Recall", "F1", "MMC", "TP", "FN",...
    "FP", "TN"];

result_table = table(analyzed_images, precision_list, recall_list,...
    F1_list, MMC_list, n_TP_List, n_FN_List, n_FP_List, n_TN_List,...
    'VariableNames', varNames);

disp(result_table)

if save
    % Save table in Excel-format
    writetable(result_table, strcat(result_save_path, ...
        'Result_Segmentation_R_Irregularis_tol=',...
        string(tolerance),'pixels.xlsx'))
end


