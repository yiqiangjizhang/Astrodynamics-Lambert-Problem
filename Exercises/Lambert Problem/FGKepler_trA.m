function [rf] = FGKepler_trA(mu, r0, r0_dot, delta_theta)

    % FGKepler_trA.m - Computes the final position vector of a spacecraft
    %
    % PROTOTYPE:
    %   [rf] = FGKepler_trA(mu, r0, r0_dot, delta_theta)
    %
    % DESCRIPTION:
    %   Computes the final position vector of a spacecraft
    %
    % INPUT:
    %   mu[1]       Gravitational parameter of the central body [km^3/s^2]
    %   r0[3]       Initial position vector of the spacecraft [km]
    %   r0_dot[3]   Initial velocity vector of the spacecraft [km/s]
    %   delta_theta[1]  Encounter angle
    %
    % OUTPUT:
    %   rf[3]       Final position vector of the spacecraft [km]
    %
    % CALLED FUNCTIONS:
    %   None
    %
    % AUTHOR:
    %   Yi Qiang Ji Zhang, 23/03/2023, MATLAB, MinETransfer.m
    %
    % -------------------------------------------------------------------------

    r0_norm = norm(r0);

    % Angular momentum
    h_vec = cross(r0, r0_dot);
    h = norm(h_vec);
    p = h ^ 2 / mu;

    sigma0 = dot(r0, r0_dot) / sqrt(mu);

    rf = zeros(3, length(delta_theta));

    % Propagate for each angle
    for i = 1:length(delta_theta)
        rf_norm = p * r0_norm / (r0_norm + (p - r0_norm) * cos(delta_theta(i)) - sqrt(p) * sigma0 * sin(delta_theta(i)));

        F = 1 - rf_norm / p * (1 - cos(delta_theta(i)));
        G = (rf_norm * r0_norm) / sqrt(mu * p) * sin(delta_theta(i));

        rf_aux = F * r0 + G * r0_dot;
        rf(1, i) = rf_aux(1);
        rf(2, i) = rf_aux(2);
        rf(3, i) = rf_aux(3);
    end

end
