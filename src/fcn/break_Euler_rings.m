function patch = break_Euler_rings(patch)

% Breaks Euler rings that are in connection to the inverstigated BP.
% Each Euler ring contains at least 3 BPs, however, if the Euler
% ring tangates the flank of the patch, one or more BPs remain undetected.

% There are x cases of interest.
% 1) The number of BPs on the ring >=3. Keep shortest geodesic path from
% investigated BP to other BPs. Remove rest of the ring.
% 2) The number of BPs on the ring = 2 (investigated BP + one other).
% Insert one acting BP, defined as the point with a geodesic distance from
% BP = D_max-1, that is further to the other BP.
% 3) The number of BPs on the ring = 1 (only the investigated BP). Insert
% two acting BPs.

% If a ring without additional BPs is a spurrious ring is impossible to
% know.

% Input: patch
% Output: patch with broken Euler ring.

[n, m] = size(patch);

BP_x = ceil(n/2);
BP_y = ceil(m/2);

% Identify all branch points
branch_points = bwmorph(patch, "branchpoints");
[bp_coords_x, bp_coords_y] = find(branch_points);
bp_coords = setdiff([bp_coords_y, bp_coords_x], [BP_y, BP_x], "rows");


% Identify the ring
test = imfill(patch, 'holes');
holes = test & ~patch;
test = test & ~branch_points;
disconnected = imfill(~test, find(holes), 4);
ring = patch & disconnected;
ring = bwareaopen(ring, 3);

% Identify the BPs on the ring.
BPs_on_ring = branch_points & ring;



if sum(BPs_on_ring(:)) == 2
    D1 = bwdistgeodesic(ring, BP_y, BP_x);
    D2 = bwdistgeodesic(ring, bp_coords(:, 1), bp_coords(:, 2));
    d_max = max(D1, [], 'all');
    [options_x, options_y] = find(D1 == (d_max - 1));
    dist_from_BP2 = [];
    for i=1:numel(options_x)
        dist_from_BP2 = [dist_from_BP2; D2(options_x(i), options_y(i))];
    end
    [~, res] = max(dist_from_BP2);
    acting_BP_x = options_x(res);
    acting_BP_y = options_y(res);
end


% Add two acting BPs
if sum(BPs_on_ring(:)) == 1
    D1 = bwdistgeodesic(ring, BP_y, BP_x);
    d_max = max(D1, [], 'all');
    half_d_max = ceil(d_max/2);
    [acting_BP_x, acting_BP_y] = find(D1 == (half_d_max+1));
end


% Broken ring
[BPs_on_ring_x, BPs_on_ring_y] = find(BPs_on_ring);
BPs_on_ring_coords = setdiff([BPs_on_ring_x, BPs_on_ring_y],...
    [BP_y, BP_x], "rows");

% Add acting BPs
try
    BPs_on_ring_coords = [BPs_on_ring_coords; [acting_BP_x, acting_BP_y]];
catch
    pass = 1;
end

% Identify part to remove
N_BP_o_R = size(BPs_on_ring_coords, 1);
part_to_remove = ring;
for i = 1:N_BP_o_R
    D1 = bwdistgeodesic(ring, BP_y, BP_x);
    D2 = bwdistgeodesic(ring, BPs_on_ring_coords(i, 2),...
        BPs_on_ring_coords(i, 1));
    D = D1 + D2;
    D = round(D * 8) / 8;
    D(isnan(D)) = inf;
    middle_path = imregionalmin(D);
    part_to_remove = part_to_remove & ~middle_path;
end

patch = patch & ~part_to_remove;

end