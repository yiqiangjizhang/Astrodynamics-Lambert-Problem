function [kep,p,E,M,dt] = car2kep(in,mu)

% car2kep.m - Convertion from Cartesian position and velocity to
%   Keplerian elements.
%
% PROTOTYPE:
%   [kep,p,E,M,dt] = car2kep(in, mu)
%   
% DESCRIPTION:
%   Convertion from Cartesian position and velocity to
%   Keplerian elements. All the units to be consistent, angles in radians.
%
% INPUT:
%	in[6]       State vector in cartesian coordinates (position [L],
%               velocity [L/T]).
%   mu          Planetary gravity constant [L^3/(M*T^2)].
%
% OUTPUT:
%   kep[1,6]    Vector of Keplerian elements: kep = [a, e, i, Om, om, th],
%               where theta is the true anomaly. a in [L],
%                   0 <= i  <= pi   [rad]
%                   0 <= Om <  2*pi [rad]
%                   0 <= om <  2*pi [rad]
%                   0 <= th <  2*pi [rad].
%   p           Parameter [L].
%   E           Eccentric anomaly, hyperbolic anomaly or parabolic anomaly
%               (for definitions see Vallado pag. 49).
%   M           Mean anomaly [rad].
%   dt          Time from the pericentre passage [T].
%
% REFERENCES:
%	D. A. Vallado, "Fundamentals of Astrodynamics and Applications, Second
%   Edition", Microcosm Press, 2001 - for computation of eccentric anomaly
%   E.
%
% CALLED FUNCTIONS:
%   crossFast
%
% ORIGINAL VERSION:
%   Massimiliano Vasile, 2002, MATLAB, ca2kep.m
%
% AUTHOR:
%   Matteo Ceriotti, 08/02/2007, MATLAB, car2kep.m
%
% PREVIOUS VERSION:
%   Matteo Ceriotti, 08/02/2007, MATLAB, cart2kep.m
%       - Header and function name in accordance with guidlines.
%
% CHANGELOG:
%   17/04/2007, REVISION: Camilla Colombo, Matteo Ceriotti
%   06/11/2007, Daniel Novak: numerical roundoff error in the calculation
%     	of theta.
%   02/12/2008, Matteo Ceriotti: Added correction for om when i == 0.
%   04/12/2008, Matteo Ceriotti: Added correction for om for planar
%     	retrograde orbits (i == pi). Changed help for bounds on Om, om, th,
%       now in [0, 2pi).
%   30/09/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%   16/10/2009, Jeannette Heiligers, Camilla Colombo:
%       - Added threshold value for eccentricity in case of circular orbit.
%       - Modified condition for determinimg is om = 2*pi-om. New condition
%           condition works for prograde/retrograde, planar/inclined case.
%       - Modified condition for determinimg is th = 2*pi-om. New condition
%           condition works for prograde/retrograde, planar/inclined case.
%           Previous condition was wrong in the circular case because
%           dot(r,v) jumps between <0 and >0.
%   19/10/2009, Matteo Ceriotti, Camilla Colombo:
%       - Added threshold value for eccentricity in case of circular orbit,
%           value determined with test_car2kep.m.
%       - Modified condition for determinimg is om = 2*pi-om. New condition
%           condition works for prograde/retrograde, planar/inclined case.
%       - Modified condition for determinimg is th = 2*pi-om. New condition
%           condition works for prograde/retrograde, planar/inclined case.
%           Previous condition was wrong in the circular case because
%           dot(r,v) jumps between <0 and >0.
%   17/11/2009, Matteo Ceriotti, Jeannette Heiligers: changed p - in the
%       output - from the 5th to the 2nd position to be consistent with
%       kep2car.
%   17/12/2009, Camilla Colombo: modified if condition to check the
%       parabolic case.
%   29/03/2010, Camilla Colombo: crossFast used.
%
% -------------------------------------------------------------------------

% Threshold on eccentricity for considering the orbit to be circular
%   Value determined comparing the relative error on state and position
%   between using the circular case and the elliptic case. With this elimit
%   the relative error on position and velocity is always less then 1e-7.
elimit = 0.00000001;

if ~isreal(in)
    error('spaceToolbox:car2kep:complexInput',...
              'Complex input');
end

r = in(1:3);
v = in(4:6);

nr = sqrt(r(1)^2+r(2)^2+r(3)^2);    % Norm of r

% Angular momentum vector: h = cross(r,v) 
h  = [r(2)*v(3)-r(3)*v(2),r(3)*v(1)-r(1)*v(3),r(1)*v(2)-r(2)*v(1)]; 
nh = sqrt(h(1)^2+h(2)^2+h(3)^2);    % Norm of h

% Inclination
i = acos(h(3)/nh);

% Line of nodes vector
if i ~= 0 && i ~= pi    % Orbit is out of xy plane
    
    % n = cross([0 0 1],h);
    % Normalisation of nv to 1: n = n/norm(n)
    %                           n = n/sqrt(n(1)^2+n(2)^2+n(3)^2);
    n = [-h(2),h(1),0]/sqrt(h(1)^2+h(2)^2);

else                    % Orbit is in xy plane
    
    % Arbitrary choice of n
    n = [1,0,0]; 
    warning('spaceToolbox:car2kep:planarOrbit','Planar orbit. Arbitrary choice of Omega = 0.');
end

% Argument of the ascending node
Om = acos(n(1));
if n(2) < 0
    Om = mod(2*pi-Om,2*pi);
end
% ---> CL: 19/10/2009, Matteo Ceriotti, Camilla Colombo: added mod(Om,2*pi)
%           to avoid numerical approximation error and hence Om=2*pi.

% Parameter
p = nh^2/mu;

% Eccentricity vector: ev = 1/mu*cross(v,h) - r/nr
ev  = 1/mu*[v(2)*h(3)-v(3)*h(2),v(3)*h(1)-v(1)*h(3),v(1)*h(2)-v(2)*h(1)] - r/nr; 
e = sqrt(ev(1)^2+ev(2)^2+ev(3)^2);    % Eccentricity (norm of eccentricity vector)

% Argument of the pericentre
if e<elimit     % Circular orbit
    
    % Arbitrary eccentricity vector
    ev = n;
    ne = 1; % ne = norm(ev)
    warning('spaceToolbox:car2kep:circularOrbit','Circular orbit. Arbitrary choice of omega = 0.');

else                    % Non circular orbit
    
    ne = e;
end         

om = acos(min(max((n(1)*ev(1) + n(2)*ev(2) + n(3)*ev(3)) / ne,-1),1)); % acos(dot(n,ev)/ne)
% In the circular case om = 0
% ---> CL: 19/10/2009, Matteo Ceriotti, Camilla Colombo: numerical roundoff
%           error in the calculation of theta to avoid argument of acos to
%           be >1 or <(-1)

if dot(h,crossFast(n,ev)) < 0 
    om = mod(2*pi-om,2*pi);
end
% ---> CL: 16/10/2009, Jeannette Heiligers, Camilla Colombo: this
%           condition works for prograde/retrograde, planar/inclined
%           case.
% ---> CL: 19/10/2009, Matteo Ceriotti, Camilla Colombo: added mod(om,2*pi)
%           to avoid numerical approximation error and hence om=2*pi.

% Semi-major axis
a = p/(1-e^2);
if a == Inf  % Parabola
    warning('spaceToolbox:car2kep:parabolicOrbit','Parabola. Semi-major axis is Inf.');
end
% ---> CL: 17/12/2009, Camilla Colombo: modified if condition to check the parabolic case.



% True anomaly: acos(dot(ev,r)/ne/nr);
th = acos(min(max((ev(1)*r(1)+ev(2)*r(2)+ev(3)*r(3))/ne/nr,-1),1));

% ---> CL: 06/11/2007, Daniel Novak: numerical roundoff error in the
%           calculation of theta to avoid argument of acos to be >1 or
%           <(-1)
% ---> CL: 16/10/2009, Jeannette Heiligers, Camilla Colombo:
%           Note: do not substitute ne with e in the previous
%           formula, because we want to have substitute [ne = norm(ev) = 1]
%           in the circular case and we want e to be the real value of the
%           eccentricity

if dot(h,crossFast(ev,r)) < 0
    % the condition dot(r,cross(h,ev)) < 0 works in the same way
    th = mod(2*pi-th,2*pi);
end
% ---> CL: 16/10/2009, Jeannette Heiligers, Camilla Colombo: this condition
%           works for prograde/retrograde, planar/inclined case
% ---> CL: 19/10/2009, Matteo Ceriotti, Camilla Colombo: added mod(th,2*pi)
%           to avoid numerical approximation error and hence th=2*pi.


kep = [a,e,i,Om,om,th];



% The following formulas have been found in Vallado chapter 2.

% ---> CL: 17/11/2009, Matteo Ceriotti, Jeannette Heiligers: changed the
%           value of nargout after providing p as the 2nd instead of the
%           5th output variable.

if nargout>2

    % E (= eccentric or hyperbolic anomaly) required as output
    if  e<elimit            % Circumference
        E = th;

    elseif e<1              % Ellipse
        sinE = sin(th)*sqrt(1-e^2)/(1+e*cos(th));
        cosE = (e+cos(th))/(1+e*cos(th));
        E = atan2(sinE,cosE);
        if E<0
            E = E+2*pi;
        end

    elseif e>1             % Hyperbola
        sinhH = sin(th)*sqrt(e^2-1)/(1+e*cos(th));
        % coshH=(e+cos(th))/(1+e*cos(th)); % Not needed
        H = asinh(sinhH);
        E = H;

    elseif e==1             % Parabola
        B = tan(th/2); % Parabolic anomaly
        E = B;

    end

    if nargout>3

        % M (= mean anomaly) required as output
        if e<elimit         % Circumference
            M=E;

        elseif e<1          % Ellipse
            M=E-e*sinE;

        elseif e>1          % Hyperbola
            M=e*sinhH-H;

        elseif e==1         % Parabola
            M=B+B^3/3;

        end

        if nargout>4

            % dt (= time from the pericentre passage) required as output
            if e<1          % Ellipse or circumference
                n=sqrt(mu/a^3);

            elseif e>1      % Hyperbola
                n=sqrt(-mu/a^3);

            elseif e==1     % Parabola
                n=2*sqrt(mu/p^3);

            end
            dt=M/n;
        end
    end
end

return