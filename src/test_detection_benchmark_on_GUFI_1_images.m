% Test


clear all
%close all
% clc


addpath(genpath('fcn')) % Add functions to path
path_img = '..\data\Ghent_University_Fungal_Images_I\Real_Images\';
path_GT = '..\data\Ghent_University_Fungal_Images_I\Hand_Labelled\';

addpath(genpath(path_img));
addpath(genpath(path_GT));

result_save_path = '..\result_tables\';
addpath(genpath(result_save_path))

% Tolerance: From Wang et al. (2019) state that they allow a tolerance of
% 2% of the image diagonal, i.e. tol = 0.02 * sqrt(2*300^2) = 
% 6*sqrt(2) = 8.4853


% Images retrived 2022-11-14

% Create image list
image_list = [];
for k=1:100
    n_digits = numel(num2str(k));
    if n_digits < 2
        image_list = [image_list; strcat("00", string(k))];
    elseif n_digits < 3
        image_list = [image_list; strcat("0", string(k))];
    else
        image_list = [image_list; string(k)];
    end
end

image_prefix = 'Crop';
image_suffux = '.png';

GT_prefix = 'MCrop';
GT_suffix = '.png';

% Vectors for storing results
precision_list = [];
recall_list = [];
F1_list = [];
MMC_list = [];
analyzed_images = [];
n_TP_List = [];
n_TN_List = [];
n_FP_List = [];
n_FN_List = [];

sigma_list = [1.5, 1.5, 1.5, 2, 1.5, 1.5, 1.5, 2, 1.5, 1.5, 1.5, 1.5,...
    2, 2, 3, 2, 2, 2, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4,...
    4, 3, 4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4,...
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,...
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4];

%% Parameters
tolerance = 8.4853; % In this paper, the tolerance is set as 2% of the length of the image diagonal.
minimum_branch_length = 5;
result_decimals = 3;
save = 1;

for i=1:100

    I = imread(strcat(image_prefix, image_list(i), image_suffux));

    % Image name for table
    analyzed_images = [analyzed_images; strcat(image_prefix,...
        image_list(i))];

    GT = imread(strcat(GT_prefix, image_list(i), GT_suffix));

    GT_bin = imbinarize(rgb2gray(GT));
    
    % Skeletonization removes redundant 4-connected pixels, 
    % 'MinBranchLength' removes small tracing errors.
    GT_bin = bwskel(GT_bin, 'MinBranchLength', 5); 


    Igray = mat2gray(I);

    [polarity, ~] = sample_polarity(Igray);

    theta = 0:15:360;

    filter_sigma = sigma_list(i);
    [sk_structure, ridgeFilt] = mycelium_detection(Igray,...
        polarity, filter_sigma, theta, minimum_branch_length);


    % Find the indices of non-zero elements
    [rowOutput, colOutput] = find(sk_structure);
    [rowGT, colGT] = find(GT_bin);

    % Initialize the cost matrix
    costMatrix = zeros(length(rowOutput), length(rowGT));

    % Calculate the Euclidean distance between all pairs of non-zero elements
    for k = 1:length(rowOutput)
        for j = 1:length(rowGT)
            dist = sqrt((rowOutput(k)-rowGT(j))^2 + ...
                (colOutput(k)-colGT(j))^2);
            costMatrix(k,j) = dist;
        end
    end


    [M,uR,uC] = matchpairs(costMatrix, tolerance/2);

    [n, m] = size(sk_structure);

    TP = size(M, 1);
    FN = size(uC, 1);
    FP = size(uR, 1);
    TN = n*m - TP - FN - FP;

    % disp(strcat("TP_test: ", string(size(M, 1))))
    % disp(strcat("FN_test: ", string(size(uR, 1))))
    % disp(strcat("FP_test: ", string(size(uC, 1))))


    precision = TP/(TP+FP);
    recall = TP/(TP+FN);
    F1 = (2*TP)/(2*TP+FP+FN);
    MMC = (TP*TN - FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)));

%     disp(strcat("Precision: ", string(precision)))
%     disp(strcat("Recall: ", string(recall)))
%     disp(strcat("F1: ", string(F1)))
%     disp(strcat("MMC: ", string(MMC)))

    % Store TP, FN, FP, and TN
    n_TP_List = [n_TP_List; TP];
    n_TN_List = [n_TN_List; TN];
    n_FP_List = [n_FP_List; FP];
    n_FN_List = [n_FN_List; FN];

    precision_list = [precision_list; precision];
    recall_list = [recall_list; recall];
    F1_list = [F1_list; F1];
    MMC_list = [MMC_list; MMC];

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


% Compute medians 
median_precision = median(precision_list);
median_recall = median(recall_list);
median_F1 = median(F1_list);
median_MMC = median(MMC_list);
median_TP = median(n_TP_List);
median_FN = median(n_FN_List);
median_FP = median(n_FP_List);
median_TN = median(n_TN_List);

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
        'Result_Detection_GUFI-1_tol=',...
        string(tolerance),'pixels.xlsx'))
end
