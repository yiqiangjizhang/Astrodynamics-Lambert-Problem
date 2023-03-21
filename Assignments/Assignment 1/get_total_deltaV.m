function [delta_V_total] = get_total_deltaV(mu_1, mu_2, r_park, r_op, v_c1, v_c2, r1_dot, rf_dot, v_dep, v_arr)

% Delta V Short
v_inf_dep = norm(r1_dot - v_dep);
v_inf_arr = norm(rf_dot - v_arr);

% Departure
v_p1 = sqrt(mu_1*(2/r_park + v_inf_dep^2/mu_1)); % Velocity at periapsis in scape orbit
delta_v1 = v_p1 - v_c1;

% Arrival
v_p2 = sqrt(mu_2*(2/r_op + v_inf_arr^2/mu_2)); % Velocity at periapsis in scape orbit
delta_v2 = v_p2 - v_c2;

% Total delta V
delta_V_total = delta_v1 + delta_v2;

end

