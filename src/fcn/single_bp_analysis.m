function [connection_point, n_arms, ridge_strengths, direction_vectors]...
    = single_bp_analysis(bp, sk_structure, ridgeFilt, patch_range,...
    theta_swich)

% Separates the intersection arms connected to a branch point (BP).
% Input:
%       bp : [1 x 2], image coordinates of a branch point.
%       sk_structure : [n x m], binary matrix repressenting the skeleton
%                      structure of the filamentous object of interest.
%       ridgeFilt : [n x m], 2D ridge filter output.
%       patch_range : [1 x 2], value specifying the size of the patch
%                     around the BP.
% Output:
%       connection_point : [1 x 2], coordinates specifying the branch out
%                          point corresponding to the BP.
%       n_arms : scalar, the number of intersection arms. 
%                n_arms = 3 -> it's a T or Y branch point.
%                n_arms = 4 -> it's a + or X branch point.
%                n_arms < 3 -> it's a false positive branch point.
%       ridge_strengths : [1 x n_arms], median ridge filter response for
%                         each arm.
%       direction_vectors : [2 x n_arms], vectors approximating the
%                           direction of each arm.

% Image size.
[n, m] = size(sk_structure);

% Define patch indices
x_range = bp(1)-patch_range(1):bp(1)+patch_range(2);
y_range = bp(2)-patch_range(1):bp(2)+patch_range(2);


% Remove invalid indices
x_range(x_range<1)=[];
x_range(x_range>m)=[];
y_range(y_range<1)=[];
y_range(y_range>n)=[];


% Extract patch
patch = sk_structure(y_range, x_range);


% Compute the corresponding coordinates of the branch points inside the
% patch.
in_patch = bp - [x_range(1) y_range(1)] + [1 1];

% Remove disconnected components.
disconnected = imfill(~patch, [in_patch(:, 2) in_patch(:, 1)], 8); 
patch = patch & disconnected;

% Euler ring case
if bweuler(patch, 8) < 1
   patch = break_Euler_rings(patch);
end

% Remore pixels surrounding BP
bp_patch = zeros(size(patch));
bp_patch(in_patch(2), in_patch(1)) = 1;
bp_patch_dil = imdilate(bp_patch, strel('square', 3));



patch_2 = patch & ~bp_patch_dil;


% Remove parts connected to other bp
other_bps = bwmorph(patch, 'branchpoints') & ~bp_patch_dil;
[x, y] = find(other_bps);

[label_patch, n_obj] = bwlabel(patch_2);


if (~isempty(x)) && (n_obj ==3)
    idx = sub2ind(size(patch), x, y);
    corrupt_arms = unique(label_patch(idx));
    for i=1:numel(corrupt_arms)
        corrupt_arm = (label_patch == corrupt_arms(i));
        corrupt_arm = corrupt_arm & ~other_bps &...
            ~remove_connection(patch, 1, [y, x]);
        [corrupt_arm_lbl, n_obj_2] = bwlabel(corrupt_arm);
        % If n_obj_2 is less than 3, do nothing.
        if n_obj_2 >= 3  
            % Identify which arm is connected.
            dist = zeros(1, n_obj_2);
            arm_lengths = zeros(1, n_obj_2);
            for j=1:n_obj_2
                part = (corrupt_arm_lbl==j);
                arm_lengths(j) = sum(part(:));
                [x, y] = find(part);
                [~, D] = dsearchn([x, y], in_patch);
                dist(j) = D;
            end
            % The part with the shortest distance to the BP is the
            % connected part.
            [~, connected_part] = min(dist);
            
            % Remove the "splitting" parts.
            indices_to_remove = setdiff(1:3, connected_part);
            parts_to_remove = (corrupt_arm_lbl == indices_to_remove(1))...
                | (corrupt_arm_lbl == indices_to_remove(2));
            label_patch(parts_to_remove) = 0;

        end
  
    end
    
end

% Extract arm coordinates
arm_coords = {n_obj};
for i=1:n_obj
    obj = (label_patch==i);
    [X, Y] = find(obj);
    arm_coords{i} = [Y, X];
end

% Extract features from all intersection arms.
n_arms = n_obj;
arm_ridge_values = zeros(1, n_arms);
arm_vectors = zeros(2, n_arms);
connection_points = zeros(2, n_arms);
closest_bps = zeros(2, n_arms);
for i=1:n_arms
    arm = arm_coords{i};
    [arm_ridge_median, arm_vector, connection_point, closest_bp] =...
        intersection_arm_features(arm, ridgeFilt, y_range, x_range,...
        in_patch);
    % Arm linreg can be poorly conditioned
    arm_ridge_values(i) = arm_ridge_median;
    arm_vectors(:, i) = arm_vector;
    connection_points(:, i) = connection_point;
    closest_bps(:, i) = closest_bp;
end

ridge_strengths = arm_ridge_values;
direction_vectors = arm_vectors;

arm_rad = atan(abs(arm_vectors(1, :))./ abs(arm_vectors(2, :)));
feature_vector = [arm_ridge_values; arm_rad];

try
    [branch_out_coords, direction_vectors] = match_branch_2(...
        patch, bp, feature_vector,...
        connection_points, x_range, y_range, direction_vectors,...
        theta_swich);
catch
    branch_out_coords = bp;
end

connection_point = branch_out_coords;
end