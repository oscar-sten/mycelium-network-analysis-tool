function tip_points = find_tip_points(sk_structure)

tips_1 = bwmorph(sk_structure, 'endpoints');
sk_structure_no_tips = sk_structure & ~tips_1;
tips_2 = bwmorph(sk_structure_no_tips, 'endpoints');
tips = tips_1 & bwareaopen((tips_1 | tips_2), 2);
[tips_x, tips_y] = find(tips);
tip_points = [tips_y, tips_x];

end