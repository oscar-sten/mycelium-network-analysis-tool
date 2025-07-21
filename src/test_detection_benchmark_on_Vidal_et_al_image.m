% Comparison between our method and the method employed in Vidal et al.
% (2019)
% https://github.com/hsueh-lab/FFT/tree/master/

% This script is very computationally demanding. 
% We recomend atleast 128 GB of RAM and Intel i9 processor or better.

clear all
%close all
% clc
addpath(genpath('fcn')) % Add functions to path
addpath(genpath('..\data\Vidal_et_al_images\'));

result_save_path = '..\result_tables\';
addpath(genpath(result_save_path))

% Images retrived 2022-11-14
image_list = ["01", "02", "03", "04" , "05", "06", "07", "08", "09", "10"];

% Add ad-hoc mask for removing ring.
% Image 01, size: 6392 x 6392


image_prefix = 'Image';
image_suffux = '.tif';

GT_prefix = 'GT_Image';
GT_suffix = '.png';

FFT_prefix = 'FFT_Image';
FFT_suffix = '.png';

precision_list_STEN = [];
precision_list_FFT = [];

recall_list_STEN = [];
recall_list_FFT = [];

F1_list_STEN = [];
F1_list_FFT = [];

MMC_list_STEN = [];
MMC_list_FFT = [];


analyzed_images = [];

n_TP_List = [];
n_TN_List = [];
n_FP_List = [];
n_FN_List = [];

%% Parameters
tolerance = 8.4853; 
result_decimals = 3;
minimum_branch_length = 3;
save = 1;

sigma_list = [2.5, 3, 3, 3.5, 3.5, 3.5, 3.5, 3.5, 3, 2.5];

disp(strcat("Tolerance: ", string(tolerance)))

for i = 1:10

    I = imread(strcat(image_prefix, image_list(i), image_suffux));
    
    disp(strcat("Processing image: ", image_list(i)));
    
    % Image name for table
    analyzed_images = [analyzed_images; strcat(image_prefix,...
        image_list(i))];

    GT = imread(strcat(GT_prefix, image_list(i), GT_suffix));

    GT_bin = imbinarize(rgb2gray(GT));

    GT_bin = bwskel(~GT_bin, 'MinBranchLength', 3);
    
    FFT_output = imread(strcat(FFT_prefix, image_list(i), FFT_suffix));

    Igray = mat2gray(I);

    [polarity, vote] = sample_polarity(Igray);

    theta = 0:15:360;
    
    filter_sigma = sigma_list(i);

    [sk_structure, ridgeFilt] = mycelium_detection(Igray,...
        polarity, filter_sigma, theta, minimum_branch_length);
    
    % Remove border
    if strcmp(image_list(i), "01") || strcmp(image_list(i), "02") ||...
            strcmp(image_list(i), "03") || strcmp(image_list(i), "09")...
            || strcmp(image_list(i), "10")
    outside = Igray<0.0001;
    outside_dil = imdilate(outside, strel('disk', 7));
    sk_structure = sk_structure & ~outside_dil;
    end


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
            %         if dist < tolerance
            %             costMatrix(i,j) = dist;
            %         end
            costMatrix(k,j) = dist;
        end
    end

    disp("Solving optimal matching problem")
     
    [M,uR,uC] = matchpairs(costMatrix, tolerance/2);

    [n, m] = size(sk_structure);

    TP = size(M, 1);
    FN = size(uC, 1);
    FP = size(uR, 1);
    TN = n*m - TP - FN - FP;


    
    precision = TP/(TP+FP);
    recall = TP/(TP+FN);
    F1 = (2*TP)/(2*TP+FP+FN);
    MMC = (TP*TN - FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)));


    precision_list_STEN = [precision_list_STEN; precision];
    recall_list_STEN = [recall_list_STEN; recall];
    F1_list_STEN = [F1_list_STEN; F1];
    MMC_list_STEN = [MMC_list_STEN; MMC];
    
    disp("Computing cost matrix for comparison")

    % Find the indices of non-zero elements
    [rowA, colA] = find(FFT_output);
    [rowB, colB] = find(GT_bin);

    % Initialize the cost matrix
    costMatrix = zeros(length(rowA), length(rowB));

    % Calculate the Euclidean distance between all pairs of non-zero elements
    for k = 1:length(rowA)
        for j = 1:length(rowB)
            dist = sqrt((rowA(k)-rowB(j))^2 + (colA(k)-colB(j))^2);
            %         if dist < tolerance
            %             costMatrix(i,j) = dist;
            %         end
            costMatrix(k,j) = dist;
        end
    end

    disp("Solving optimal matching problem for comparison")
    [M,uR,uC] = matchpairs(costMatrix, tolerance/2);

    [n, m] = size(sk_structure);

    TP = size(M, 1);
    FN = size(uC, 1);
    FP = size(uR, 1);
    TN = n*m - TP - FN - FP;
    precision = TP/(TP+FP);
    recall = TP/(TP+FN);
    F1 = (2*TP)/(2*TP+FP+FN);
    MMC = (TP*TN - FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)));


    precision_list_FFT = [precision_list_FFT; precision];
    recall_list_FFT = [recall_list_FFT; recall];
    F1_list_FFT = [F1_list_FFT; F1];
    MMC_list_FFT = [MMC_list_FFT; MMC];

    disp("Done")
end


%% Summarize results in a table

% Compute averages
average_precision = mean(precision_list_STEN);
average_recall = mean(recall_list_STEN);
average_F1 = mean(F1_list_STEN);
average_MMC = mean(MMC_list_STEN);


% Add averages to tables
analyzed_images = [analyzed_images; "Average"];
precision_list_STEN = [precision_list_STEN; average_precision];
recall_list_STEN = [recall_list_STEN; average_recall];
F1_list_STEN = [F1_list_STEN; average_F1];
MMC_list_STEN = [MMC_list_STEN; average_MMC];


% Round to decired number of decimals
precision_list_STEN = round(precision_list_STEN, result_decimals);
recall_list_STEN = round(recall_list_STEN, result_decimals);
F1_list_STEN = round(F1_list_STEN, result_decimals);
MMC_list_STEN = round(MMC_list_STEN, result_decimals);

%% FFT

% Compute averages
average_precision = mean(precision_list_FFT);
average_recall = mean(recall_list_FFT);
average_F1 = mean(F1_list_FFT);
average_MMC = mean(MMC_list_FFT);

% Add averages to tables
precision_list_FFT = [precision_list_FFT; average_precision];
recall_list_FFT = [recall_list_FFT; average_recall];
F1_list_FFT = [F1_list_FFT; average_F1];
MMC_list_FFT = [MMC_list_FFT; average_MMC];

% Round to decired number of decimals
precision_list_FFT = round(precision_list_FFT, result_decimals);
recall_list_FFT = round(recall_list_FFT, result_decimals);
F1_list_FFT = round(F1_list_FFT, result_decimals);
MMC_list_FFT = round(MMC_list_FFT, result_decimals);


%% Differences
precision_diff = precision_list_STEN - precision_list_FFT;
recall_diff = recall_list_STEN - recall_list_FFT;
F1_diff = F1_list_STEN - F1_list_FFT;
MMC_diff = MMC_list_STEN - MMC_list_FFT;


% Create table
varNames = ["Image", "Precision STEN", "Precision FFT", "Precision Diff",...
    "Recall STEN", "Recall FFT", "Recall Diff" "F1 STEN", "F1 FFT",...
    "F1 Diff", "MMC STEN", "MMC FFT", 'MMC Diff'];

result_table = table(analyzed_images, precision_list_STEN,...
    precision_list_FFT, precision_diff, recall_list_STEN, recall_list_FFT,...
    recall_diff, F1_list_STEN, F1_list_FFT, F1_diff,...
    MMC_list_STEN, MMC_list_FFT, MMC_diff,...
    'VariableNames', varNames);

disp(result_table)


%% Save
if save
    % Save table in Excel-format
    writetable(result_table, strcat(result_save_path, ...
        'Result_Segmentation_FFT_images_A_Oligospora_tol=', ...
        string(tolerance), 'pixels.xlsx'))
end



