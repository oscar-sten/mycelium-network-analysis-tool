function [arm_ridge_median, arm_vector, connection_point, closest_bp] =...
    intersection_arm_features(arm, ridgeFilt, y_range, x_range,...
    in_patch)

% Input: 
%        arm: the arm coordinates in image patch coordinates.
%        ridgeFilt: ridge filter output, a matrix.
%        y_range, x_range: specififying the relation between image patch
%        coordinates and global image coordinates.
%        in_patch: BP in image patch coordinates.
% Output: 
%        arm_ridge_median: median ridge filer response in arm coordinates.
%        arm_vector: vector repressenting arm in a cartestian coordinate
%                    system with the origin in the center of the image
%                    patch.
%        connection_point: connection point, i.e., the point in the arm
%                          closest to the BP, in image patch coordinates.
%        closest_bp: closest BP to connection point in image patch
%                    coordinates.

% Transform to global image coordinates.
arm_X_global = arm(:, 2) + y_range(1)-1;
arm_Y_global = arm(:, 1) + x_range(1)-1;

% Comupute median ridge value of arm coordinates.
arm_length = size(arm_X_global, 1);
arm_ridge_resp = zeros(1, arm_length);
for j=1:arm_length
    arm_ridge_resp(j) = ridgeFilt(arm_X_global(j),...
        arm_Y_global(j));
end
arm_ridge_median = median(arm_ridge_resp);

% Identify the connection point, i.e., the point in the arm closest to the
% BP.
[K,D] = dsearchn(arm, in_patch);
[~, ind] = min(D);
ind_closest = K(ind);
connection_point = arm(ind_closest, :);
closest_bp = in_patch(ind, :);

% Affine transfromation for transforming from image to cartesian
% coordinates. 
carth_origin = [ceil(numel(x_range)/2) ceil(numel(y_range)/2)];
A = [1 0;
    0 -1];

b = [-carth_origin(1);
    carth_origin(2)];

% Transform arm point to cartesian coordinates.
arm = A*arm' + b;

% Perform linear regression.
line = polyfit(arm(1, :), arm(2, :), 1);

% Compute points estimating the arm vector line.
min_x = min(arm(1, :));
max_x = max(arm(1, :));
x_range_2 = (min_x:max_x);
line_points_y = -((x_range_2.*(-line(1)))-line(2));
arm_line_estimate = [x_range_2' line_points_y'];

% Transform connection point to cartesian coordinates.
connection_point_T = A*connection_point' + b;

% Extract the point in the line furthest from the connection point
[~, ind_furthest] = max(vecnorm(arm_line_estimate...
    - arm(:, ind_closest)', 2, 2));

% Compure the arm vector as the difference between the furthest point on
% the line estimating the arm and the connection point (in cartesian 
% coordinates).
arm_vector = arm_line_estimate(ind_furthest, :) - connection_point_T';

% Transpose arm vector
arm_vector = arm_vector';


end