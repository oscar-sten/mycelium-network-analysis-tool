function removed_connection = remove_connection(patch, n_points, in_patch)

% If an intersectionpoint if connected through a right angle, i.e.,
    %  0 1 0            0 1 0
    %  0 1 1    ->      0 0 1
    %  0 0 0            0 0 0
    % the diagonal connection needs to be broken for bwlabel to work
    % correctly, thus they have to first be removed, then put back with
    % correct label.

% Input:
%       patch : [X x Y], patch of binary skeleton matrix.
%       n_points : number of points of interest.
%       in_patch : [n_points x 2], in-patch coordinates of points of
%                  interest.
% 
% Output:
%       remove_connection : [X x Y], binary matrix with the indices to be
%                           removed as 1.
    
removed_connection = zeros(size(patch));
for i=1:n_points
    test = in_patch(i, :);
    try
    if (patch(test(2)+1, test(1)) ||...
            patch(test(2)-1, test(1))) && ...
            (patch(test(2), test(1)+1) || ...
            patch(test(2), test(1)-1))
        removed_connection(test(2)+1, test(1)) = patch(test(2)+1,...
            test(1));
        removed_connection(test(2)-1, test(1)) = patch(test(2)-1,...
            test(1));
        removed_connection(test(2), test(1)+1) = patch(test(2),...
            test(1)+1);
        removed_connection(test(2), test(1)-1) = patch(test(2),...
            test(1)-1);
    end
    catch
        % if the index is at the border, just ignore
    end
end


end