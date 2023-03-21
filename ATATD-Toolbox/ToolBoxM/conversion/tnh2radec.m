function x_radec = tnh2radec(x_tnh)

% tnh2radec.m - Vector reference frame transformation:
%   Tangential-normal-h to Right Ascension and DEClination reference frame.
%
% PROTOTYPE:
%   x_radec = tnh2radec(x_tnh)
%
% DESCRIPTION:
%   Transformation from tangential-normal-h to right ascension and
%   declination reference frame.
%   Tangential-normal-h (tnh) reference frame: {t,n,h}
%       t-axis: direction tangent to the motion.
%       h-axis: direction of the angular momentum.
%       n-axis: direction in the orbit plane, normal to t (inward).
%   Right ascension and declination -spherical equatorial- (radec)
%   reference frame: {r,alpha,delta}
%       r       modulus of the vector [L].
%       alpha	in-plane right ascention angle, counted from the
%               tangential direction to the projection of the vector on the
%               orbital plane [rad].
%       delta 	out-of-plane declination angle from the projection of the
%               vector on the orbital plane up to the vector itself [rad].
%   Note: other definition use the {t,h,n} reference frame, where n is
%       outward!
%
% INPUT:
%	x_tnh[3]        Vector to be transformed, expressed in {t,n,h}.
%
% OUTPUT:
%	x_radec[3,1]    Vector transformed into {r,alpha,delta} (column).
%
% CALLED FUNCTIONS:
%   (none)
%
% AUTHOR:
%   Camilla Colombo, 23/10/2007, MATLAB, tnh2radec.m
%
% PREVIOUS VERSION:
%   Camilla Colombo, 23/10/2007, MATLAB, tnh_radecT.m
%       - Header and function name in accordance with guidlines.
%   
% CHANGELOG:
%   12/12/2007, REVISION: Matteo Ceriotti
%   11/02/2008, Matteo Ceriotti: Help improved.
%   26/03/2009, Camilla Colombo: in case x_radec(1) = 0, x_radec(2) and 
%       x_radec(3) are set to 0 as default.
%   30/09/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------

x_radec = [0;0;0];

x_radec(1) = (x_tnh(1)^2+x_tnh(2)^2+x_tnh(3)^2)^0.5;
if x_radec(1) == 0
    return
end

x_radec(2) = atan2(x_tnh(2),x_tnh(1));
x_radec(3) = asin(x_tnh(3)/x_radec(1));

return