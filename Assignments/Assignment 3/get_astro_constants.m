function [mu_Sun, mu_Earth, R_Sun, R_Earth] = get_astro_constants()

% Mu from planets [km3/s2]
mu_Sun = getAstroConstants('Sun', 'mu');
mu_Earth = getAstroConstants('Earth', 'mu');

% Radius of planets [km]
R_Sun = getAstroConstants('Sun', 'Radius');
R_Earth = getAstroConstants('Earth', 'Radius');

end

