function [mu_Sun, mu_Earth, mu_Mars, R_Sun, R_Earth, R_Mars] = get_astro_constants()

% Mu from planets [km3/s2]
mu_Sun = getAstroConstants('Sun', 'mu');
mu_Earth = getAstroConstants('Earth', 'mu');
mu_Mars = getAstroConstants('Mars', 'mu');

% Radius of planets [km]
R_Sun = getAstroConstants('Sun', 'Radius');
R_Earth = getAstroConstants('Earth', 'Radius');
R_Mars = getAstroConstants('Mars', 'Radius');

end

