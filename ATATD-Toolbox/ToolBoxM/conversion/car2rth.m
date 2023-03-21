function x_rth = car2rth(x_car,s_car)

% car2rth.m - Vector reference frame transformation:
%   Cartesian to radial-trasversal-h reference frame.
%
% PROTOTYPE:
%   x_rth = car2rth(x_car, s_car)
%
% DESCRIPTION:
%   Transformation from Cartesian to radial-trasversal-h reference frame.
%   Cartesian (car) reference frame: {x,y,z}
%       inertial reference frame.
%   Radial-transversal-h (rth) reference frame: {r,t,h}
%       r-axis: direction of the orbit radius.
%       h-axis: direction of the angular momentum.
%       t-axis: in the orbit plane, completes the reference frame. If the
%               orbit is circular it is in the direction of the velocity
%               vector.
%
% INPUT:
%	x_car[3]    Vector to be transformed, expressed in {x,y,z}.
%	s_car[6]    State vector (position [L], velocity [L/T]) of the orbiting
%               body, expressed in {x,y,z}.
%
% OUTPUT:
%	x_rth[3,1]  Vector transformed into {r,t,h}.
%
% EXAMPLE:
%	Given a spacecraft in orbit:
%       - we have the thrust vector in {x,y,z};
%       - we want the thrust vector in {r,t,h}.
%   In this case:
%       x_car = Thrust vector in {x,y,z};
%       s_car = [position, velocity] of the spacecraft in {x,y,z};
%       x_rth = Thrust vector, transformed into {r,t,h}.
%
% CALLED FUNCTIONS:
%   crossFast
%
% AUTHOR:
%   Camilla Colombo, 03/03/2006, MATLAB, car2rth.m
%
% PREVIOUS VERSION:
%   Camilla Colombo, 03/03/2006, MATLAB, car_rthT.m
%       - Header and function name in accordance with guidlines.
%   
% CHANGELOG:
%   10/01/2007, REVISION: Matteo Ceriotti
%   24/01/2008, Matteo Ceriotti, Nicolas Croisard: Help improved.
%   30/09/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%   29/03/2010, Camilla Colombo: crossFast used.
%
% -------------------------------------------------------------------------

x_car = x_car(:);
s = s_car(:);
r = s(1:3);
v = s(4:6);
r_ = r/norm(r);
h = crossFast(r,v);
h_ = h/norm(h);
t_ = crossFast(h_,r_);

A = [r_ t_ h_];

x_rth = A'*x_car;

return