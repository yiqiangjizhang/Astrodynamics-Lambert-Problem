function x_radec = car2radec(x_car,s_car)

% car2radec.m - Vector reference frame transformation:
%	Cartesian to Right Ascension and DEClination reference frame.
%
% PROTOTYPE:
%   x_radec = car2radec(x_car, s_car)
%
% DESCRIPTION:
%   Transformation from Cartesian to right ascension and declination
%   reference frame.
%   Cartesian (car) reference frame: {x,y,z}
%       inertial reference frame
%   Right ascension and declination -spherical equatorial- (radec)
%   reference frame: {r,alpha,delta}
%       r       modulus of the vector [L].
%       alpha	in-plane right ascention angle, counted from the
%               tangential direction to the projection of the vector on the
%               orbital plane [rad].
%       delta 	out-of-plane declination angle from the projection of the
%               vector on the orbital plane up to the vector itself [rad].
%
% INPUT:
%	x_car[3]        Vector to be transformed, expressed in {x,y,z}.
%	s_car[3]        State vector (position [L], velocity [L/T]) of the
%                   orbiting body, expressed in {x,y,z}.
%
% OUTPUT:
%	x_radec[3,1]    Vector transformed into {r,alpha,delta}.
%
% EXAMPLE:
%   Given a spacecraft in orbit:
%       - we have the thrust vector in {x,y,z};
%       - we want the thrust vector in {r,alpha,delta}.
%   In this case:
%       x_car = Thrust vector in {x,y,z};
%       s_car = [position, velocity] of the spacecraft in {x,y,z};
%       x_radec = Thrust vector, transformed into {r,alpha,delta}.
%
% CALLED FUNCTIONS:
%   car2tnh, tnh2radec.
%
% AUTHOR:
%   Camilla Colombo, 23/10/2007, MATLAB, car2radec.m
%
% PREVIOUS VERSION:
%   Camilla Colombo, 23/10/2007, MATLAB, car_radecT.m
%       - Header and function name in accordance with guidlines.
%   
% CHANGELOG:
%   12/12/2007, REVISION: Matteo Ceriotti
%   11/02/2008, Matteo Ceriotti: Help improved.
%   30/09/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------

x_tnh = car2tnh(x_car,s_car);
x_radec = tnh2radec(x_tnh);

return