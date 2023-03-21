function x_tnh = car2tnh(x_car,s_car)

% car2tnh.m - Vector reference frame transformation:
%   Cartesian to tangential-normal-h reference frame.
%
% PROTOTYPE:
%   x_tnh = car2tnh(x_car, s_car)
%
% DESCRIPTION:
%   Transformation from Cartesian to tangential-normal-h reference frame.
%   Cartesian (car) reference frame: {x,y,z}
%       inertial reference frame.
%   Tangential-normal-h (tnh) reference frame: {t,n,h}
%       t-axis: direction tangent to the motion.
%       h-axis: direction of the angular momentum.
%       n-axis: direction in the orbit plane, normal to t (inward).
%   Note: other definition use the {t,h,n} reference frame, where n is
%       outward!
%
% INPUT:
%	x_car[3]    Vector to be transfromed, expressed in {x,y,z}.
%	s_ca[6]     State vector (position [L], velocity [L/T]) of the orbiting
%               body, expressed in {x,y,z}.
%
% OUTPUT:
%	x_tnh[3,1]  Vector transformed into {t,n,h}.
%
% EXAMPLE:
%   Given a spacecraft in orbit:
%       - we have the thrust vector in {x,y,z};
%       - we want the thrust vector in {r,n,h}.
%   In this case:
%       x_car = Thrust vector in {x,y,z};
%       s_car = [position, velocity] of the spacecraft in {x,y,z};
%       x_tnh = Thrust vector, transformed into {t,n,h}.
%
% CALLED FUNCTIONS:
%   crossFast
%
% AUTHOR:
%   Camilla Colombo, 03/03/2006, MATLAB, car2tnh.m
%
% PREVIOUS VERSION:
%   Camilla Colombo, 03/03/2006, MATLAB, car_tnhT.m
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

x_car = x_car(:);
s = s_car(:);
r = s(1:3);
v = s(4:6);
t_ = v/norm(v);
h = crossFast(r,v);
h_ = h/norm(h);
n_ = crossFast(h_,t_);

A = [t_ n_ h_];

x_tnh = A'*x_car;

return