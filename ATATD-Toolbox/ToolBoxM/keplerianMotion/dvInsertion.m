function dv = dvInsertion(s,t,vf,rp,e)

% dvInsertion.m - Computes the dv to insert in an orbit around a planet.
%   
% PROTOTYPE:
%   dv = dvInsertion(s, t, vf, rp, e)
%
% DESCRIPTION:
%   It computes the dv to insert in an orbit around a planet. It assumes a
%   patched conics framework, and performs the manoeuvre at the pericentre
%   of the incoming hyperbola to achieve the target orbit.
%
% INPUT:
%   s[1]	Planet (1 <= s <= 9).
%   t[1]    Time of the insertion (final time of the trajectory)
%           [d, MJD2000].
%   vf[3]   Absolute cartesian velocity at the planet arrival [km/s].
%   rp[1]   Radius of the pericentre of the orbit to achieve around the
%           planet [km].
%   e[1]    Eccentricity of the orbit to achieve.
%
% OUTPUT:
%   dv[1]   Delta-v needed to insert into the defined orbit [km/s].
%
% CALLED FUNCTIONS:
%   astroConstants, ephSS_car
%
% AUTHOR:
%   Matteo Ceriotti, 19/05/2007, MATLAB, dvInsertion.m
%
% PREVIOUS VERSION:
%   Matteo Ceriotti, 19/05/2007, MATLAB, dv_insertion.m
%       - Header and function name in accordance with guidlines.
%
% CHANGELOG:
%   21/05/2007, REVISION: Camilla Colombo
%   30/12/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------

vf = vf(:);

mu = astroConstants(10+s);
[xp2,vp2] = ephSS_car(s,t);
vp2 = vp2';

vinf = norm(vf-vp2);        % Velocity relative to the planet (at infinite on the hyperbola)
vp1 = sqrt(vinf^2+2*mu/rp); % Velocity at the pericentre of the hyperbola
vp2 = sqrt(2*mu/rp-mu/rp*(1-e)); % Velocity at the pericentre of the target orbit
dv = abs(vp1-vp2);          % Delta-v at the pericentre

return