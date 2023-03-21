function [kep,mass,M] = ephNEO(time,id)

% ephNEO.m - Ephemerides of Near Earth Objects
%
% PROTOTYPE:
%	[kep, mass, M] = ephNEO(time, id)
%
% DESCRIPTION:
%   This function returns the orbital parameters, the mass, and the mean
%   anomaly of some NEOs. Each NEO is identified by an id.
%
% INPUT:
%   time[1]     MJD2000 [d].
%   id[1]       NEO identifier. It is a number starting from 12 (because
%               the identifiers from 1 to 11 are reserved for the Solar
%               System).
%
% OUTPUT:
%   kep[1,6]    Keplerian parameters. It is a 6 entry vector:
%               [a e i Om om wom] where:
%                   a is the semimajor axis [km];
%                   e is the eccentricity;
%                   i is the inclination [rad];
%                   Om is the anomaly of the ascending node [rad];
%                   om is the anomaly of the pericentre [rad];
%                   wom is the true anomaly (from the pericentre) [rad].
%   mass        Mass of the NEO [kg]. It can be read from the database, or,
%               if not available, estimated by an approximate equation.
%   M           Mean anomaly at time [rad].
%
% CALLED FUNCTIONS:
%   astroConstants
%
% AUTHOR:
%   Paolo De Pascale,  November 2004, MATLAB, ephNEO.m
%
% PREVIOUS VERSION:
%   Paolo De Pascale,  November 2004, MATLAB, NeoEphemeris.m
%       - Header and function name in accordance with guidlines.
%
% CHANGELOG:
%   20/06/2006, Christie Maddock, Camilla Colombo, Pau Sanchez
%   20/06/2006, Matteo Ceriotti: Ccleaned up code, added references to
%       astro_constants.
% 	07/07/2006, Pau Sanchez: Added data on more NEOs.
%   25/08/2006, Pau Sanchez: Added data on more NEOs.
% 	12/02/2007, Pau Sanchez: Modified mass Itokawa from 6.01e10 kg to
%       3.51e10 kg - updated from papers.
% 	16/05/2007, Camilla Colombo: added new PHO form Paolo and Max database
%                      Note for Paolo and Max: id = 57 to 72 of your
%                      code is here id 68 to 83.
% 	07/06/2008, Camilla Colombo: added 2002AA29.
%   04/10/2010, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%   18/02/2019, Joan-Pau Sanchez: Remove all old asteroids, and uploaded
%   matrix for Advanced Topics of Astrodynamics and Trajectory Design.
%
% -------------------------------------------------------------------------
persistent neo

if isempty(neo)
%     load('PHAEphemeris_ATATD_2018_2019.mat');
%     neo=[PHANEOsATATD20182019 zeros(length(PHANEOsATATD20182019(:,8)),1)];
    load('allNEOEphemeris_ATATD_2018_2019.mat');
    neo=[allNEOEphemeris_ATATD_2018_2019 zeros(length(allNEOEphemeris_ATATD_2018_2019(:,8)),1)];
end

id=round(id)-11;


% to UPDATE use this layout:
% %   name, class
% neo(id,1)=semimajor axis [AU]     neo(id,2)=eccentricity
% neo(id,3)=inclination [deg]     neo(id,4)=asc. node/raan [deg]
% neo(id,5)=arg. perigee [deg]     neo(id,6)=mean anomoly, M at time given in neo(id,8) [deg]
% neo(id,7)=abs magnitude (i.e. intrinsic brightness)     neo(id,8)=time at which Mo is given [MJD]   neo(id,9)=mass [kg]


%=====================================================

AU  =   astroConstants(2);         % km
mu_sun  =   astroConstants(4);         % km^3/s^2 Sun gravitational constant

a  =neo(id,1)*AU;
e  =neo(id,2);
di =neo(id,3);
omm=neo(id,4);
om =neo(id,5);
m  =neo(id,6);
n  =sqrt(mu_sun/a^3);

t0 =m*pi/180/n;

timediff=neo(id,8)-51544.5; % Convert to MJD2000
t=(time-timediff)*86400+t0;


p=2*pi*sqrt(a^3/mu_sun);
np=floor(t/p);
t=t-p*np;
phi=n*t;
M  = phi;
if M>pi
    M = M-2*pi;
end

% Ciclo di Newton: devo mandare a zero la funzione
%f=(M-E+esinE)...parto da E=M...deltaE=-fdot^-1*f
for i=1:5
    ddf=(e*cos(phi)-1);
    phi=phi-t*n/ddf+(phi-e*sin(phi))/ddf;
end

wom=2*atan(sqrt((1+e)/(1-e))*tan(phi*0.5)); % In radians
kep=[a e di*pi/180 omm*pi/180 om*pi/180 wom];

if neo(id,9) ~= 0
    mass=neo(id,9);
else
    %   polynomial fit to the H-d diagram which gives the diameter of the NEO as a function of its magnitude.
    %   polynomial of fifth order, Note: NOT VERY ACCURATE!
    d=-2.522e-2*neo(id,7)^5+3.2961*neo(id,7)^4-1.7249e2*neo(id,7)^3+4.5231e3*neo(id,7)^2-5.9509e4*neo(id,7)+3.1479e5;
    %   estimated mass of the NEO computed considering a density of 2 kg/dm^3
    mass=4*pi/3*(0.5*d)^3*2*1000;
end

return
