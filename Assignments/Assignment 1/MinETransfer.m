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
%   a_min[1] 
%   e_min[1] 
%   r[1,3]   Cartesian position of the body (Sun-centered for all bodies,
%            Earth-centered for the Moon) [km].
%   v[1,3]   Cartesian velocity [km/s].
%   t_m[1]   Transfer Method. Value is either "+1" for Short Way or "-1" for
%            Long Way
%
% OUTPUT:
%   a_min[1] Minimum semi-major axis [km]
%   e_min[1] Minimum eccentricity [km]
%   delta_t_min[1] Minimum Transfer time [sec]
%   r[1,3]   Cartesian position of the body (Sun-centered for all bodies,
%            Earth-centered for the Moon).
%   v[1,3]   Cartesian velocity 
%   r1_dot[1,3]   Cartesian velocity of the spacecraft
%   delta_theta[1] Encounter angle
%

c = norm(r2 - r1);
r1_norm = norm(r1);
r2_norm = norm(r2);

a_min = 0.25*(r1_norm + r2_norm + c);

cos_theta = (dot(r1,r2))/(r1_norm*r2_norm);
sin_theta = t_m * sqrt(1 - cos_theta^2);
delta_theta = acos(cos_theta);


p_min = (r1_norm*r2_norm)/c * (1 - cos_theta);
e_min = sqrt(1 - p_min/a_min);
beta_e = 2 * asin(sqrt((2*a_min - c)/(2*a_min)));

F = 1 - r2_norm/p_min*(1 - cos_theta);
G = (r2_norm*r1_norm)/(sqrt(mu*p_min)) * sin_theta;
        
r1_dot = G^(-1) * (r2 - F*r1);
delta_t_min = sqrt(a_min^3/mu) * (pi - t_m*(beta_e - sin(beta_e)));

end