function x_tnh = radec2tnh(x_radec)

% radec2tnh.m - Vector reference frame transformation.
%   Right Ascension and DEClination to tangentential-normal-h reference
%   frame.
%
% PROTOTYPE:
%	x_tnh = radec2tnh(x_radec)
%
% DESCRIPTION:
%   Transformation from right ascension and declination to
%   tangential-normal-h reference frame.
%   Right ascension and declination -spherical equatorial- (radec)
%   reference frame: {r,alpha,delta}
%       r       modulus of the vector [L].
%       alpha	in-plane right ascention angle, counted from the
%               tangential direction to the projection of the vector on the
%               orbital plane [rad].
%       delta 	out-of-plane declination angle from the projection of the
%               vector on the orbital plane up to the vector itself [rad].
%   Tangential-normal-h (tnh) reference frame: {t,n,h}
%       t-axis: direction tangent to the motion.
%       h-axis: direction of the angular momentum.
%       n-axis: direction in the orbit plane, normal to t (inward).
%   Note: other definition use the {t,h,n} reference frame, where n is
%       outward!
%
% INPUT:
%	x_radec[3]  Vector to be transformed, expressed in {r,alpha,delta}.
%
% OUTPUT:
% 	x_tnh[3,1]  Vector transformed into {t,n,h} (column).
%
% CALLED FUNCTIONS:
%   (none)
%
% AUTHOR:
%   Camilla Colombo, 07/12/2007, MATLAB, radec2tnh.m
%
% PREVIOUS VERSION:
%   Camilla Colombo, 07/12/2007, MATLAB, radec_tnhT.m
%       - Header and function name in accordance with guidlines.
%   
% CHANGELOG:
%   12/12/2007, REVISION: Matteo Ceriotti
%   11/02/2008, Matteo Ceriotti: Help improved.
%   30/09/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------

x_tnh = [0;0;0];

x_tnh(1) = x_radec(1) * cos(x_radec(3))*cos(x_radec(2));
x_tnh(2) = x_radec(1) * cos(x_radec(3))*sin(x_radec(2));
x_tnh(3) = x_radec(1) * sin(x_radec(3));

return