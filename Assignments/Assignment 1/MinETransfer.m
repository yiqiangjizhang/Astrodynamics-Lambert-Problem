function [a_min, e_min, delta_t_min, r1_dot, delta_theta] = MinETransfer(mu, r1, r2, t_m)

    % MinETransfer.m - Calculate the Minimum Energy Transfer
    %
    % PROTOTYPE:
    %   [a_min, e_min, delta_t_min] = MinETransfer(r_dep_Earth, r_arr_Mars, t_m);
    %
    % DESCRIPTION:
    %   Computes the Minimum Energy Transfer using general approach to
    %   Lambert's problem
    %
    % INPUT:
    %   mu[1]    Gravitational parameter of the central body [km^3/s^2]
    %   r1[1,3]  Cartesian position of the departure body [km]
    %   r2[1,3]  Cartesian position of the arrival body [km]
    %   t_m[1]   Transfer Method. Value is either "+1" for Short Way or "-1" for
    %            Long Way
    %
    % OUTPUT:
    %   a_min[1]        Minimum semi-major axis [km]
    %   e_min[1]        Minimum eccentricity [unitless]
    %   delta_t_min[1]  Minimum Transfer time [sec]
    %   r1_dot[1,3]     Cartesian velocity of the spacecraft at departure [km/s]
    %   delta_theta[1]  Encounter angle
    %
    % CALLED FUNCTIONS:
    %   None
    %
    % AUTHOR:
    %   Yi Qiang Ji, 23/03/2023, MATLAB, MinETransfer.m
    %
    % -------------------------------------------------------------------------

    % Distance between the two bodies
    c = norm(r2 - r1);
    r1_norm = norm(r1);
    r2_norm = norm(r2);

    % Minimum semi-major axis
    a_min = 0.25 * (r1_norm + r2_norm + c);

    cos_theta = (dot(r1, r2)) / (r1_norm * r2_norm);
    sin_theta = t_m * sqrt(1 - cos_theta ^ 2);
    delta_theta = acos(cos_theta);

    p_min = (r1_norm * r2_norm) / c * (1 - cos_theta);
    e_min = sqrt(1 - p_min / a_min);
    beta_e = 2 * asin(sqrt((2 * a_min - c) / (2 * a_min)));

    % F and G constants
    F = 1 - r2_norm / p_min * (1 - cos_theta);
    G = (r2_norm * r1_norm) / (sqrt(mu * p_min)) * sin_theta;

    % Velocity at departure
    r1_dot = G ^ (-1) * (r2 - F * r1);
    delta_t_min = sqrt(a_min ^ 3 / mu) * (pi - t_m * (beta_e - sin(beta_e)));

end
