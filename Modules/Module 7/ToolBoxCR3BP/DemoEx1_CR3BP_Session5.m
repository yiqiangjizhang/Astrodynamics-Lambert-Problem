%% DemoEx1_CR3BP_Session4.m
% NOTES:
% This exercise allows comparing the shape of an elliptic orbit in both the
% inertial reference frame and the rotating reference frame.
% 1. Run each of the sections:
%    Click on each section in order and click Ctrl+Enter to run
%% PART 0. Set up
%--------------------------------------------------------------------------
clear all; clc; close all
%--------------------------------------------------------------------------
% 1.1 Constants and definitions for Earth Moon System
massEarth   = getAstroConstants('Earth','mass');
muEarth   = getAstroConstants('Earth','mu');
massMoon = getAstroConstants('Moon','mass');
mu = massMoon/(massMoon+massEarth);
G = astroConstants(1);
EarthMoonDistance=379700; % km
RAD = pi/180;
RE = getAstroConstants('Earth','Radius');
%--------------------------------------------------------------------------
%% PART 1. Define orbit and plot in inertial reference frame 
%--------------------------------------------------------------------------
% Plot the frame
UnitDistance=EarthMoonDistance;
run PlotEarthMoon_system_Inertialframe.m
%--------------------------------------------------------------------------
% Define a Keplerian Orbit
rp = 0.0840; % periapsis distance (in Earth_Moon distance units)
ra = 0.75;   % apoapsos distance (in Earth_Moon distance units)

% rp = 0.4; % periapsis distance (in Earth_Moon distance units)
% ra = 0.86;   % apoapsos distance (in Earth_Moon distance units)
% 
% rp = 0.2; % periapsis distance (in Earth_Moon distance units)
% ra = 0.7615;   % apoapsos distance (in Earth_Moon distance units)

sma = (rp+ra)/2;       % semimajor axis
ecc = (ra-rp)/(rp+ra); % eccentricity 
inc = 0;  % inclination (rad)
OM0 = 0;  % RAAN (rad)
om  = pi; % Argument of the periapsis (rad)
trA = 0;  % true anomaly (rad)
%--------------------------------------------------------------------------
% Propagate Keplerian Orbit
% x.1. Compute cartesian inertial state vector
KepAsteroid = [sma ecc inc OM0 om trA];
PAst=2*pi*sqrt(sma^3/(1-mu));
SV = kep2car(KepAsteroid,1-mu);
% x.2 Propagate orbit
TSPAN = [0 4*PAst];
OPTIONS = odeset('RelTol',3e-10,'AbsTol',1e-12);  % lower accuracy
MODEL = @dyn_2BP;
[t,Xinert]     = ode45(MODEL,TSPAN,SV,OPTIONS,1-mu);
%--------------------------------------------------------------------------
% Plot Keplerian Orbit
  plot3(Xinert(:,1),Xinert(:,2),Xinert(:,3))
%% PART 2. Plot same orbit in rotating reference frame 
%--------------------------------------------------------------------------
% Plot the frame
run PlotEarthMoon_system_frame.m
%--------------------------------------------------------------------------
% x.2. Normalize units to CR3BP
X=SV(1:3);
V=SV(4:6);
% x.3. Is a rotation necessary? nop because at t=0, both Earth and asteroid
% are at trA=0
% x.4. Account for the mu translation of the center
X=[X(1)-mu X(2) X(3)]; % moved to the baricentre
V=[V(1) V(2)-mu V(3)];
% x.5. Rotating Reference frame
X=[X(1) X(2) X(3)];
V=[V(1)+X(2) V(2)-X(1) V(3)];
%--------------------------------------------------------------------------
% Propagate in CR3BP 3D
OPTIONS = odeset('RelTol',3e-10,'AbsTol',1e-12);

TintR = [0 4*PAst];
[SV_forward,t_forward] = trajGet3BP_3D([X V],0,TintR(2),mu,OPTIONS);
%--------------------------------------------------------------------------
plot3(SV_forward(:,1),SV_forward(:,2),SV_forward(:,3))
