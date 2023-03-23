function [delta_V_total] = get_total_deltaV(mu_1, mu_2, r_park, r_op, v_c1, v_c2, r1_dot, rf_dot, v_dep, v_arr)

    % get_total_deltaV.m - Computes the total delta V of a transfer
    %
    % PROTOTYPE:
    %   [delta_V_total] = get_total_deltaV(mu_1, mu_2, r_park, r_op, v_c1, v_c2, r1_dot, rf_dot, v_dep, v_arr)
    %
    % DESCRIPTION:
    %   Computes the total delta V of a transfer
    %
    % INPUT:
    %   mu_1[1]        Gravitational parameter of the departure body [km^3/s^2]
    %   mu_2[1]        Gravitational parameter of the arrival body [km^3/s^2]
    %   r_park[1]      Periapsis radius of the parking orbit [km]
    %   r_op[1]        Periapsis radius of the operational orbit [km]
    %   v_c1[1]        Velocity of the departure body [km/s]
    %   v_c2[1]        Velocity of the arrival body [km/s]
    %   r1_dot[3]      Velocity vector of the spacecraft at the departure body [km/s]
    %   rf_dot[3]      Velocity vector of the spacecraft at the arrival body [km/s]
    %   v_dep[3]       Velocity vector of the spacecraft at the departure body [km/s]
    %   v_arr[3]       Velocity vector of the spacecraft at the arrival body [km/s]

    %
    % OUTPUT:
    %   delta_V_total[1]    Total delta V of the transfer [km/s]
    %
    % CALLED FUNCTIONS:
    %   None
    %
    % AUTHOR:
    %   Yi Qiang Ji, 23/03/2023, MATLAB, get_total_deltaV.m
    %
    % -------------------------------------------------------------------------

    % Delta V Short
    v_inf_dep = norm(r1_dot - v_dep);
    v_inf_arr = norm(rf_dot - v_arr);

    % Departure
    v_p1 = sqrt(mu_1 * (2 / r_park + v_inf_dep ^ 2 / mu_1)); % Velocity at periapsis in scape orbit
    delta_v1 = v_p1 - v_c1;

    % Arrival
    v_p2 = sqrt(mu_2 * (2 / r_op + v_inf_arr ^ 2 / mu_2)); % Velocity at periapsis in scape orbit
    delta_v2 = v_p2 - v_c2;

    % Total delta V
    delta_V_total = delta_v1 + delta_v2;

end
