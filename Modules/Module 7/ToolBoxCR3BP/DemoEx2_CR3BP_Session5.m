%% DemoEx2_CR3BP_Session4.m
% NOTES:
clear all; close all; clc
%% Part 1 Compute Collinear equilibrium positions 
% 1.1. Define the system and compute the mass parameter mu
    massEarth   = getAstroConstants('Earth','mass');
    massMoon = getAstroConstants('Moon','mass');
    mu = massMoon/(massMoon+massEarth);
EarthMoonDistance=379700; % km
  
%--------------------------------------------------------------------------
% Solve explicit condition for: L1
    fun=@(x)x-(1-mu)/(x+mu)^2+mu/(x-1+mu)^2; % Definition 
    x0 = 0.5; %first Guess
    xL1 = fzero(fun,x0);     
%--------------------------------------------------------------------------  
%% Part 2 - plot trajectorories for equilibrium positions
%--------------------------------------------------------------------------
% 2.1 - Propagate state vector for 1 full period of the Moon

    SV0=[xL1 0 0 0 0 0]; % state vector of equilibrium position

    OPTIONS = odeset('RelTol',3e-10,'AbsTol',1e-12);
    MODEL = @pcr3bp_3D;
    TintR = [0 2*pi];
    [times,SV] = ode45(MODEL,TintR,SV0,OPTIONS,mu,[]);

%--------------------------------------------------------------------------
% 2.2 Plot the resulting trajectory in rotating reference frame
UnitDistance=EarthMoonDistance;
run PlotEarthMoon_system_frame.m
plot3(SV(:,1),SV(:,2),SV(:,3),'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
    'Marker','pentagram')
%--------------------------------------------------------------------------
% 2.3 Use function synodic2car to convert each point in the trajectory to
% inertial reference frame.

theta_0=0; % Initial angular position of the Moon
X_c=zeros(length(times),3);
V_c=zeros(length(times),3);

for it=1: length(times)
    [X_c(it,:), V_c(it,:)] = ...
        synodic2car(SV(it,1:3), SV(it,1:3), times(it), mu, theta_0);
end

%--------------------------------------------------------------------------
% 2.4. Plot the resulting trajectory in inertial reference frame
UnitDistance=EarthMoonDistance;
run PlotEarthMoon_system_Inertialframe.m
plot3(X_c(:,1),X_c(:,2),X_c(:,3))