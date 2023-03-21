function x_car = rth2car(x_rth,s_car)

% rth2car.m - Vector reference frame transformation.
%   Radial-trasversal-h to Cartesian reference frame.
%
% PROTOTYPE:
%	x_car = rth2car(x_rth, s_car)
%
% DESCRIPTION:
%   Transformation from radial-trasversal-h to Cartesian reference frame.
%   Radial-transversal-h (rth) reference frame: {r,t,h}
%       r-axis: direction of the orbit radius.
%       h-axis: direction of the angular momentum.
%       t-axis: in the orbit plane, completes the reference frame. If the
%               orbit is circular it is in the direction of the velocity
%               vector.
%   Cartesian (car) reference frame: {x,y,z}
%       inertial reference frame.
%
% INPUT:
%	x_rth[3]    Vector to be transformed, expressed in {r,t,h}.
%  	s_car[6]    State vector (position [L], velocity [L/T]) of the orbiting
%               body, expressed in {x,y,z}.
%
% OUTPUT:
%  	x_car[3,1]  Vector transformed into {x,y,z}.
%
% EXAMPLE:
%   Given a spacecraft in orbit:
%       - we have the thrust vector in {r,t,h};
%       - we want the thrust vector in {x,y,z}.
%   In this case:
%       x_rth = Thrust vector in {r,t,h};
%       s_car = [position, velocity] of the spacecraft in {x,y,z};
%       x_car = Thrust vector, transformed in {x,y,z}.
%
% CALLED FUNCTIONS:
%   crossFast
%
% AUTHOR:
%   Camilla Colombo, 03/03/2006, MATLAB, rth2car.m
%
% PREVIOUS VERSION:
%   Camilla Colombo, 03/03/2006, MATLAB, rth_carT.m
%       - Header and function name in accordance with guidlines.
%   
% CHANGELOG:
%   10/01/2007, REVISION: Matteo Ceriotti
%   11/02/2008, Matteo Ceriotti: Help improved.
%   30/09/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%   29/03/2010, Camilla Colombo: crossFast used.
%
% -------------------------------------------------------------------------

x_rth = x_rth(:);
s = s_car(:);
r = s(1:3);
v = s(4:6);
r_ = r/norm(r);
h = crossFast(r,v);
h_ = h/norm(h);
t_ = crossFast(h_,r_);

A = [r_ t_ h_];

x_car = A*x_rth;

return