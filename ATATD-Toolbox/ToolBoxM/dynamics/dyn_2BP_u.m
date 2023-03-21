function sp = dyn_2BP_u(t,s,u,mu)

% dyn_2BP_u.m - Equation of the Dynamics of the Two-Body Problem with
%   acceleration control vector u in Cartesian reference frame.
%
% PROTOTYPE:
%   sp = dyn_2BP_u(t, s, u, mu)
%
% DESCRIPTION:
%   Computes the derivative of the state vector (position and velocity) for
%   the two-body problem dynamics with control acceleration, in Cartesian
%   reference frame.
%
% INPUT:
%	t[1]    Time [T].
%	s[6]    State vector (position, velocity) in Cartesian reference frame
%           [L L/T].
%  	u[3]    Control acceleration vector in Cartesian coordinates [L/T^2].
%   mu[1]   Planetary constant [L^3/(M*T^2)].
%           mu = G*mp with  mp = mass of the attracting body
%                           G  = universal gravity constant.
%
% OUTPUT:
% 	sd      Derivative of the state vector [L/T L/T^2]
%
% CALLED FUNCTIONS:
%   (none)
%
% AUTHOR:
%   Camilla Colombo, 07/03/2005, MATLAB, dyn_2BP_u.m
%
% PREVIOUS VERSION:
%   Camilla Colombo, 07/03/2005, MATLAB, d2b_eq_u.m
%       - Header and function name in accordance with guidlines.
%
% CHANGELOG:
%   10/01/2007, REVISION: Matteo Ceriotti
%   04/01/2010, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------

u = u(:);

sp = [       s(4:6)       ;...
      -mu/norm(s(1:3))^3*s(1:3) + u];

return