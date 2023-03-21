function x_car = radec2car(x_radec,s_car)

% radec2car.m - Vector reference frame transformation:
%	Right Ascension and DEClination to Cartesian reference frame.
%
% PROTOTYPE:
%   x_car = radec2car(x_radec, s_car)
%
% DESCRIPTION:
%   Transformation from right ascension and declination to Cartesian
%   reference frame.
%   Right ascension and declination -spherical equatorial- (radec)
%   reference frame: {r,alpha,delta}
%       r       modulus of the vector [L].
%       alpha	in-plane right ascention angle, counted from the
%               tangential direction to the projection of the vector on the
%               orbital plane [rad].
%       delta 	out-of-plane declination angle from the projection of the
%               vector on the orbital plane up to the vector itself [rad].
%   Cartesian (car) reference frame: {x,y,z}
%       inertial reference frame
%
% INPUT:
%	x_radec[3]  Vector to be transformed, expressed in {r,alpha,delta}.
% 	s_car[6]    State vector (position [L], velocity [L/T]) of the orbiting
%               body, expressed in {x,y,z}.
%
% OUTPUT:
% 	x_car[3,1]  Vector transformed into {x,y,z}.
%
% EXAMPLE:
%   Given a spacecraft in orbit:
%       - we have the thrust vector in {r,alpha,delta};
%       - we want the thrust vector in {x,y,z}.
%   In this case:
%       x_radec = Thrust vector in {r,alpha,delta};
%       s_car = [position, velocity] of the spacecraft in {x,y,z};
%       x_car = Thrust vector, transformed in {x,y,z}.
%
% CALLED FUNCTIONS:
%   radec2tnh, tnh2car.
%
% AUTHOR:
%   Camilla Colombo, 07/12/2007, MATLAB, radec2car.m
%
% PREVIOUS VERSION:
%   Camilla Colombo, 07/12/2007, MATLAB, radec_carT.m
%       - Header and function name in accordance with guidlines.
%   
% CHANGELOG:
%   12/12/2007, REVISION: Matteo Ceriotti
%   11/02/2008, Matteo Ceriotti: Help improved.
%   30/09/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------

x_tnh = radec2tnh(x_radec);
x_car = tnh2car(x_tnh,s_car);

return