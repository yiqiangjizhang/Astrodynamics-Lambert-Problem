function rho = lrom(kep_rob,doe)

% lrom.m - Linearized Relative Orbit Motion for general elliptic orbits.
%
% PROTOTYPE:
%   rho = lrom(kep_rob, doe)
%
% DESCRIPTION:
%   It calculates the relative position coordinates in terms of the orbit
%   element differences. They are expressed in a radial-transverse-h
%   {r,theta,h} reference frame.
%
% INPUT:
% 	kep_rob[6]	Orbital parameters of the chief orbit [a e i Om om f]
%             	f is the true anomaly at which the relative position is
%             	computed [rad]. Units are: a [L], Om, om , i, f [rad].
% 	doe[6]      Orbit element differences between the chief orbit and the
%               perturbed orbit [da de di dOm dom dM], where for example
%               da = a_perturbed - a_chief. dM is the difference in mean
%               anomaly at the point at which the the relative position is
%             	computed [rad]. Units are: da [L], dOm, dom , di, dM [rad].
%
% OUTPUT:
% 	rho[3,1]    Relative position vector expressed in {r,theta,h} reference
%               frame:
%                   r-axis: direction of the orbit radius.
%                   h-axis: direction of angular momentum.
%                   theta-axis: transverse direction in the orbit plane,
%                           completes the reference frame. If the
%                           orbit is circular it is in the direction of the
%                           velocity vector.
%
% CALLED FUNCTIONS:
%   (none)
%
% REFERENCES:
%	- Schaub H., Junkins J. L., "Analytical Mechanics of Space Systems",
%     AIAA Educationa Series, 2003, pp. 619 - Relative Motion Equations.
%
% AUTHOR:
%   Camilla Colombo, 07/06/2006, MATLAB, lrom.m
%
% PREVIOUS VERSION:
%   Camilla Colombo, 07/06/2006, MATLAB, lrom.m
%       - Header and function name in accordance with guidlines.
%
% CHANGELOG:
%   19/05/2007, REVISION: Matteo Ceriotti
%   29/03/2010, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------

a  = kep_rob(1);     da  = doe(1);
e  = kep_rob(2);     de  = doe(2);
i  = kep_rob(3);     di  = doe(3);
Om = kep_rob(4);     dOm = doe(4);
om = kep_rob(5);     dom = doe(5);
f  = kep_rob(6);      
                     dM  = doe(6);

eta = sqrt(1-e^2);
r = a*(1-e^2)/(1+e*cos(f));
theta = om+f;

% Coordinate relative to the chief orbit: x = x(f)
%                                         y = y(f)
%                                         z = z(f)

x = r/a*da + a*e*sin(f)/eta*dM - a*cos(f)*de;
y = r/eta^3*(1+e*cos(f))^2*dM + r*dom + r*sin(f)/eta^2*(2+e*cos(f))*de + r*cos(i)*dOm;
z = r*(sin(theta)*di - cos(theta)*sin(i)*dOm);
rho = [x y z]';

%--------------------------------------------------------------------------

% Linearized non-dimensional relative orbit motion

% fu = atan(e*dM/(-eta*de));
% fv = fu-pi/2;   % atan(eta*de/(e*dM))
% thetaw = atan(di/(-sin(i)*dOm));
% du = sqrt(e^2*dM^2/eta^2 + de^2);
% dw = sqrt(di^2 + sin(i)^2*dOm^2);

% u = da/a - e*de/(2*eta^2) + du/eta^2*(cos(f-fu)+e/2*cos(2*f-fu));
% v = ((1+e^2/2)*dM/eta^3 + dom + cos(i)*dOm) - du/eta^2*(2*sin(f-fu)+e/2*sin(2*f-fu));
% w = dw*cos(theta-thetaw);
% rho = [u v w]';

%--------------------------------------------------------------------------

return