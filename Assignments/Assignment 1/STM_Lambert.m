function [STM] = STM_Lambert(mu, rf, delta_E, F, G, r0, r0_dot, delta_t)

    % STM_Lambert.m - Computes the State Transition Matrix for Lambert's problem
    %
    % PROTOTYPE:
    %   [STM] = STM_Lambert(mu, rf, delta_E, F, G, r0, r0_dot, delta_t)
    %
    % DESCRIPTION:
    %   Computes the State Transition Matrix for Lambert's problem
    %
    % INPUT:
    %   mu[1]           Gravitational parameter [km^3/s^2]
    %   rf[3]           Cartesian position of the spacecraft at arrival [km]
    %   delta_E[1]      Eccentric anomaly [rad]
    %   F[1]            Kepler constant F
    %   G[1]            Kepler constant G
    %   r0[3]           Cartesian position of the spacecraft at departure [km]
    %   r0_dot[1,3]     Cartesian velocity of the spacecraft at departure [km/s]
    %   delta_t[1]      Transfer time [sec]

    %
    % OUTPUT:
    %   STM[3,3]        State Transition Matrix
    %
    % CALLED FUNCTIONS:
    %   None
    %
    % AUTHOR:
    %   Yi Qiang Ji, 23/03/2023, MATLAB, STM_Lambert.m
    %
    % -------------------------------------------------------------------------

    % Normalize vectors
    r0_norm = norm(r0);
    r0_dot_norm = norm(r0_dot);
    rf_norm = norm(rf);

    % Semi-major axis
    a = mu / (2 * mu / r0_norm - r0_dot_norm ^ 2);

    F_dot =- (sqrt(mu * a)) / (rf_norm * r0_norm) * sin(delta_E);
    G_dot = 1 - a / rf_norm * (1 - cos(delta_E));

    rf_dot = F_dot * r0 + G_dot * r0_dot;

    C = a * sqrt(a ^ 3 / mu) * (3 * sin(delta_E) - (2 + cos(delta_E)) * delta_E) - ...
        a * delta_t * (1 - cos(delta_E));

    delta_r = rf - r0;
    delta_v = rf_dot - r0_dot;

    STM = r0_norm / mu * (1 - F) * (delta_r' * r0_dot - delta_v' * r0) + ...
        C / mu * rf_dot' * r0_dot + G * eye(3);
end
