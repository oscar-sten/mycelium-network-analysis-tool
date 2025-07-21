function is_crossing = crossing_criterion(dist, dist_limit, phi_bp_1...
    , phi_bp_2, theta_T)

min_angle = max(abs(phi_bp_1-90), abs(phi_bp_2-90));

% Boolean function
is_crossing = (min_angle > theta_T) || (dist < dist_limit);

end