function v1 = eulerAxisAngle(v,n,theta)

% eulerAxisAngle.m - Euler axis and angle rotation.
%
% PROTOTYPE:
%   v1 = eulerAxisAngle(v, n, theta)
%
% DESCRIPTION:
%   Rotates a vector about an axis of a given angle (counterclockwise
%   according to the right-hand rule).
%   Note: If you want to rotate the coordinate system of a given angle
%   (counterclockwise according to the right-hand rule), use the minus in
%   front of the angle (i.e., -theta). (See reference).
%
% INPUT:
%   v[3]        Vector to be rotated.
%   n[3]        Axis of rotation.
%   theta       Angle of rotation [rad].
%
% OUTPUT:
%   v1[3,1]     Rotated vector.
%
% REFERENCES:
%	Schaub and Junkins, "Analytical Mechanics of Space Systems", AIAA
%	Education Series, 2003, pp. 90.
%   see http://mathworld.wolfram.com/RotationMatrix.html for Matrix, vector
%   rotation.
%
% CALLED FUNCTIONS:
%   (none)
%
% AUTHOR:
%   Matteo Ceriotti, 11/01/2007, MATLAB, eulerAxisAngle.m
%
% PREVIOUS VERSION:
%   Matteo Ceriotti, 11/01/2007, MATLAB, euler_axis_angle.m
%       - Header and function name in accordance with guidlines.
%   
% CHANGELOG:
%   14/05/2007, REVISION: Camilla Colombo
%   30/09/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%   28/03/2011, Camilla Colombo, Matteo Ceriotti: Corrected code and
%       description accordingly. In order to rotate a vector about an axis
%       of a given angle (counterclockwise according to the right-hand rule)
%       the matrix R must be transposed.
%       Note: eulerAxisAngle now rotates a vector about an axis of a given
%       angle (counterclockwise according to the right-hand rule) or it
%       rotates the coordinate system of a given (clockwise)
%       angle.
%
% -------------------------------------------------------------------------

v = v(:);
n = n/norm(n);
R = [cos(theta)+(1-cos(theta))*n(1)^2, (1-cos(theta))*n(1)*n(2)+sin(theta)*n(3), (1-cos(theta))*n(1)*n(3)-sin(theta)*n(2);...
        (1-cos(theta))*n(1)*n(2)-sin(theta)*n(3), cos(theta)+(1-cos(theta))*n(2)^2, (1-cos(theta))*n(2)*n(3)+sin(theta)*n(1);...
        (1-cos(theta))*n(1)*n(3)+sin(theta)*n(2), (1-cos(theta))*n(2)*n(3)-sin(theta)*n(1), cos(theta)+(1-cos(theta))*n(3)^2];

% --> CL 28/03/2011, Camilla Colombo, Matteo Ceriotti. Matrix R transposed.
v1 = R'*v;

return