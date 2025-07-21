% Manual length measurement

clear all
%close all
% clc
addpath(genpath('fcn')) % Add functions to path
addpath(genpath('..\data\lbl_length_Cardini_et_al_images\'));

%% Parameters
plotting = 0;
disconnectedComponentsRemovalFlag = 1;
filter_sigma = 4;
minimum_branch_length = 20;


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

    if strcmp(test_image_list(i), 'h')
        mm_per_pixel = 0.0011;
    else
        [Igray, pixels_scale_bar, mm_scale_bar] = remove_scale_bar_2(...
            Igray);

        mm_per_pixel = mm_scale_bar/pixels_scale_bar;
    end

    %% Compute auto
    [polarity, vote] = sample_polarity(Igray);
    theta = [0:15:360];
    [sk_structure, ridgeFilt] = mycelum_segmentation(Igray,...
        polarity, filter_sigma, theta, minimum_branch_length);
    
    %% Remove disconnected components
    if disconnectedComponentsRemovalFlag
        [sk_structure] = remove_disconnected_components(sk_structure);
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
                    pt_1 = data_3(k, :);
                    pt_2 = data_3(k+1, :);
                    dist = norm(pt_1-pt_2);
                    pixel_dist = pixel_dist + dist;
                end

            end
        end
    end

    if plotting
        [Y, X] = find(sk_structure);
        scatter(X, Y, 'Marker', '.', 'MarkerEdgeColor', 'red')
    end

    length_mm = mm_per_pixel*pixel_dist;
    length_auto_mm = compute_length(sk_structure, mm_per_pixel);
    if plotting
        title(strcat("Length manual: ", string(round(length_mm, 2)),...
            "mm, Length auto: ", string(round(length_auto_mm, 2)), "mm"))
    end

    manual_length_list = [manual_length_list; length_mm];
    auto_length_list = [auto_length_list; length_auto_mm];

end

%% 
figure()
plot(manual_length_list, auto_length_list, 'bo', 'LineWidth', 1)
hold on
% Using GLM here is pretty much useless, but the result is the same as for
% a least squares or chi 2 regression
glm = fitglm(manual_length_list, auto_length_list);
coeff = glm.Coefficients.Estimate;
r2 = glm.Rsquared.Ordinary;
x_range = linspace(0, 80);
plot(x_range, coeff(1) + coeff(2)*x_range, 'red', 'LineWidth', 1)
plot(x_range, x_range, 'LineStyle', '--', 'Color', 'black', 'LineWidth', 1)

xlabel('Manual measurements [mm]', 'Interpreter','Latex')
ylabel('Automatic measurements [mm]', 'Interpreter','Latex')
axis([0 80 0 80])
legend('Results', strcat("y = ",string(round(coeff(1), 2)), ...
    " + ",string(round(coeff(2), 2)), 'x', ", R^2 = ",...
    string(round(r2, 2))), 'Line of agreement', 'Location', 'northwest')

% Lin's Concordance Correlation Coefficient
CCC = f_CCC(...
    [manual_length_list, auto_length_list], 0.05);
CCC_hat = CCC{1}.est;

title(strcat("Estimated length of R. irregulris in-vivo, $\hat{\rho_c}$ = ",...
    string(round(CCC_hat, 2))), 'Interpreter','Latex')

lower_conf = CCC{1}.confInterval(1);

if lower_conf < 0.9
    assessment = "poor";
elseif lower_conf > 0.90 && lower_conf < 0.95
    assessment = "moderate";
elseif lower_conf > 0.95 && lower_conf < 0.99
    assessment = "substantial";
elseif lower_conf > 0.99 && lower_conf < 1
    assessment = "almost perfect";
else
    error("Invalid concordance value")
end

disp(strcat("Lower 95% confidence interval is: ",...
    string(round(lower_conf, 2)),...
    ", which according to McBride (2005) can motivate the assessment that the concordance is: ",...
    assessment))






