% Create GT Skeletons

clear all
%close all
% clc
addpath(genpath('fcn')) % Add functions to path

addpath(genpath('..\data\lbl_length_Cardini_et_al_images\'));

save_path = '..\data\ground_truth_skeletons_Cardini_et_al_images\';
addpath(genpath(save_path))

plotting = 0;
save = 1;

test_image_list = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'l',...
    'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v'];

suffix_img = '.jpg';
suffix_annotation = '_lines_lbl.mat';

n_annotated_images = numel(test_image_list);

manual_length_list = [];
auto_length_list = [];

for i = 1:n_annotated_images
    I = imread(strcat(test_image_list(i), suffix_img));

    Igray = rgb2gray(I);

    GT_skeleton = zeros(size(Igray));

    if strcmp(test_image_list(i), 'h')
        mm_per_pixel = 0.0011;
    else
        [Igray, pixels_scale_bar, mm_scale_bar] = remove_scale_bar_2(...
            Igray);

        mm_per_pixel = mm_scale_bar/pixels_scale_bar;
    end

%% Manual

    I_lbl = load(strcat(test_image_list(i), suffix_annotation));
    label_data = I_lbl.gTruth.LabelData;
    VarNames = label_data.Properties.VariableNames;
    all_points = [];
    pixel_dist = 0;
    if plotting
        figure()
        imshow(I)
        hold on
    end
    
    for v=VarNames
        data = label_data.(v{1});
        for d = 1:numel(data)
            data_2 = data{d};
            for j = 1:numel(data_2)
                data_3 = data_2{j};
                all_points = [all_points; data_3];
                if plotting
                    plot(data_3(:, 1), data_3(:, 2), 'b')
                end
                [n_pts, ~] = size(data_3);
                for k = 1:n_pts-1
                    % Extract data points
                    pt_1 = data_3(k, :);
                    pt_2 = data_3(k+1, :);

                    % Add skeleton line
                    [n, m] = size(GT_skeleton);
                    nPoints = max(abs(pt_1(1) - pt_2(1)), ...
                        abs(pt_1(2) - pt_2(2)))+1;
                    rIndex = round(linspace(pt_1(2), ...
                        pt_2(2), nPoints));  % Row indices
                    cIndex = round(linspace(pt_1(1),...
                        pt_2(1), nPoints));  % Column indices
                    if rIndex(1)<1
                        rIndex(1) = [];
                        cIndex(1) = [];
                    end
                    if cIndex(1)<1
                        rIndex(1) = [];
                        cIndex(1) = [];
                    end
                    if rIndex(end)>n
                        rIndex(end) = [];
                        cIndex(end) = [];
                    end
                    
                    if cIndex(end)>m
                        cIndex(end) = [];
                        rIndex(end) = [];
                    end


                    index = sub2ind([n, m], rIndex, cIndex);
                    GT_skeleton(index) = 1;
                    
                    % Compute distance
                    dist = norm(pt_1-pt_2);
                    pixel_dist = pixel_dist + dist;
                end

            end
        end
    end


if save
    imwrite(GT_skeleton, strcat(save_path, "GT_skeleton_",...
        test_image_list(i), '.png'))
end






end
