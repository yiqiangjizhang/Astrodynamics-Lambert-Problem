function [rf, rf_dot, delta_E_sol, F, G] = FGKepler_dt(mu, r0, r0_dot, delta_t)

    % FGKepler_dt.m - Function to compute the final position and velocity of a spacecraft
    %
    % PROTOTYPE:
    %   [rf, rf_dot, delta_E_sol, F, G] = FGKepler_dt(mu, r0, r0_dot, delta_t)
    %
    % DESCRIPTION:
    %   This function computes the final position and velocity of a spacecraft and the F and G constants
    %
    % INPUT:
    %   mu[1]           Gravitational parameter [km^3/s^2]
    %   r0[3]           Initial position vector [km]
    %   r0_dot[3]       Initial velocity vector [km/s]
    %   delta_t[1]      Transfer time [sec]
    %
    % OUTPUT:
    %   rf[3]           Final position vector [km]
    %   rf_dot[3]       Final velocity vector [km/s]
    %   delta_E_sol[1]  Eccentric anomaly [rad]
    %   F[1]            F constant
    %   G[1]            G constant

    %
    % CALLED FUNCTIONS:
    %   None
    %
    % AUTHOR:
    %   Yi Qiang Ji Zhang, 23/03/2023, MATLAB, FGKepler_dt.m
    %
    % -------------------------------------------------------------------------

    % Compute the initial position and velocity norm
    r0_norm = norm(r0);
    r0_dot_norm = norm(r0_dot);

    % Compute the semi-major axis
    a = mu / (2 * mu / r0_norm - r0_dot_norm ^ 2);

    n = sqrt(mu / a ^ 3);

    delta_M = n * delta_t;

    sigma_0 = (dot(r0, r0_dot)) / sqrt(mu);

    % Solve for the eccentric anomaly
    fun = @(delta_E) delta_E - (1 - r0_norm / a) * sin(delta_E) - sigma_0 / sqrt(a) * (cos(delta_E) - 1) - delta_M; % function
    x0 = delta_M; % initial point
    delta_E_sol = fzero(fun, x0);

    % Check
    fun(delta_E_sol);

    F = 1 - a / r0_norm * (1 - cos(delta_E_sol));
    G = delta_t + sqrt(a ^ 3 / mu) * (sin(delta_E_sol) - delta_E_sol);

    rf = F * r0 + G * r0_dot;
    rf_norm = norm(rf);

    G_dot = 1 - a / rf_norm * (1 - cos(delta_E_sol));

    rf_dot = 1 / G * (G_dot * rf - r0);

end
