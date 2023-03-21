function x_tnh = rth2tnh(x_rth,a,e,f,mu)

% rth2tnh.m - Vector reference frame transformation.
%   Radial-transversal-h to tangential-normal-h reference frame.
%
% PROTOTYPE:
%	x_tnh = rth2tnh(x_rth, a, e, f, mu)
%
% DESCRIPTION:
%   Transformation from radial-trasversal-h to tangential-normal-h
%   reference frame.
%   Radial-transversal-h (rth) reference frame: {r,t,h}
%       r-axis: direction of the orbit radius.
%       h-axis: direction of the angular momentum.
%       t-axis: in the orbit plane, completes the reference frame. If the
%               orbit is circular it is in the direction of the velocity
%               vector.
%   Tangential-normal-h (tnh) reference frame: {t,n,h}
%       t-axis: direction tangent to the motion.
%       h-axis: direction of the angular momentum.
%       n-axis: direction in the orbit plane, normal to t (inward).
%   Note: other definition use the {t,h,n} reference frame, where n is
%       outward!
%
% INPUT:
%	x_rth[3]    Vector to be transformed, expressed in {r,t,h}.
%	a           Semi-major axis [L].
%  	e           Eccentricity.
% 	f           True anomaly from the pericentre [rad].
%   mu          Planetary gravity constant [L^3/(M*T^2)].
%
% OUTPUT:
%  	x_tnh[3,1]  Vector transformed into {t,n,h}.
%
% EXAMPLE:
%   Given a spacecraft in orbit:
%       - we have the thrust vector in {r,t,h};
%       - we want the thrust vector in {t,n,h}.
%   In this case:
%       x_rth = Thrust vector in {r,t,h};
%       a,e,f = Orbital parameters of the spacecraft;
%       mu = Gravitational constant of the attractor;
%       x_tnh = Thrust vector, transformed in {t,n,h}.
%
% CALLED FUNCTIONS:
%   (none)
%
% AUTHOR:
%   Camilla Colombo, 10/03/2006, MATLAB, rth2tnh.m
%
% PREVIOUS VERSION:
%   Camilla Colombo, 10/03/2006, MATLAB, rth_tnhT.m
%       - Header and function name in accordance with guidlines.
%   
% CHANGELOG:
%   07/02/2006, Camilla Colombo: removed v from inputs.
%   10/01/2007, REVISION: Matteo Ceriotti
%   11/02/2008, Matteo Ceriotti: Help improved.
%   30/09/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------

x_rth = x_rth(:);
p = a*(1-e^2);
n = sqrt(mu/a^3);
h = n*a^2*sqrt(1-e^2);
r = p/(1+e*cos(f));
v = sqrt(2*mu/r - mu/a);

sinb = h*e/(p*v)*sin(f);
cosb = h/(p*v)*(1+e*cos(f));

x_tnh = [sinb -cosb 0; cosb sinb 0; 0 0 1]'*x_rth;

return;