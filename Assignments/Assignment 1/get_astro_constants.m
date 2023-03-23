function [mu_Sun, mu_Earth, mu_Mars, R_Sun, R_Earth, R_Mars] = get_astro_constants()

    % get_astro_constants.m - Astronomical constants of celestial bodies.
    %
    % PROTOTYPE:
    %	mu_Sun, mu_Earth, mu_Mars, R_Sun, R_Earth, R_Mars] = get_astro_constants()
    %
    % DESCRIPTION:
    %   It calculates the gravitational parameter and radius of the Sun, Earth and Mars.
    %
    % INPUT:
    %   None
    %
    % OUTPUT:
    %   mu_Sun[1]       Gravitational parameter of the Sun [km3/s2]
    %   mu_Earth[1]     Gravitational parameter of the Earth [km3/s2]
    %   mu_Mars[1]      Gravitational parameter of Mars [km3/s2]
    %   R_Sun[1]        Radius of the Sun [km]
    %   R_Earth[1]      Radius of the Earth [km]
    %   R_Mars[1]       Radius of Mars [km]
    %
    % CALLED FUNCTIONS:
    %   None
    %
    % AUTHOR:
    %   Yi Qiang Ji, 23/03/2023, MATLAB, get_astro_constants.m
    %
    % -------------------------------------------------------------------------

    % Mu from planets [km3/s2]
    mu_Sun = getAstroConstants('Sun', 'mu');
    mu_Earth = getAstroConstants('Earth', 'mu');
    mu_Mars = getAstroConstants('Mars', 'mu');

    % Radius of planets [km]
    R_Sun = getAstroConstants('Sun', 'Radius');
    R_Earth = getAstroConstants('Earth', 'Radius');
    R_Mars = getAstroConstants('Mars', 'Radius');

end
