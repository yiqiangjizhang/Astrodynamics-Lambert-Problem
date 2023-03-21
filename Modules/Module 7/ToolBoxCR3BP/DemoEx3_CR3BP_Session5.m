%% DemoEx4_CR3BP_Session4.m
% NOTES: Compute the zero velocity curves for a trajectory with Jacobi
% Constant C=L1.
clear all; close all; clc
%% Part 1 Compute Collinear equilibrium positions 
% 1.1. Define the system and compute the mass parameter mu
    massEarth   = getAstroConstants('Earth','mass');
    massMoon = getAstroConstants('Moon','mass');
    mu = massMoon/(massMoon+massEarth);
EarthMoonDistance=379700; % km
    G = astroConstants(1);
%--------------------------------------------------------------------------
% Solve explicit condition for: L1
    fun=@(x)x-(1-mu)/(x+mu)^2+mu/(x-1+mu)^2; % Definition 
    x0 = 0.5; %first Guess
     xL1 = fzero(fun,x0);     
%--------------------------------------------------------------------------  
%% Part 2 - Compute Jacobi Constant of L1 equilibrium
%--------------------------------------------------------------------------
% normalization units
%   unitLength = EarthMoonDistance;
%     muSystem = G*(massMoon+massEarth); 
%     P        = 2*pi*sqrt(unitLength^3/muSystem);
%     normTime = P/2/pi;
%     normVel  = unitLength/normTime; 
%--------------------------------------------------------------------------
    SV0=[xL1 0 0 0 0 0]; % state vector of equilibrium position

    x  = SV0(1);
    y  = SV0(2);
    z  = SV0(3);
    vx = SV0(4);
    vy = SV0(5);
    vz = SV0(6);
    xEarth = -mu;
    xMoon  = 1-mu;
    rho1=sqrt((x-xEarth)^2+y^2+z^2);
    rho2=sqrt((x-xMoon)^2+y^2+z^2);    
    
    v=sqrt(vx^2+vy^2+vz^2);
    C1=x^2+y^2+2*(1-mu)/rho1+2*mu/rho2-v^2;
%--------------------------------------------------------------------------    
%% Ploting zero velocity curves in 3D
%--------------------------------------------------------------------------  
delta=-1.5:0.01:1.5;
[X, Y, Z]=meshgrid(delta, delta, delta);

rho1=sqrt((X+mu).^2 + Y.^2 + Z.^2);
rho2=sqrt((X-(1-mu)).^2 + Y.^2 + Z.^2);

V2=(X.^2+Y.^2)+2*(1-mu)./rho1+2*mu./rho2-C1;
%--------------------------------------------------------------------------
% 2.2 Plot the resulting trajectory in rotating reference frame
UnitDistance=EarthMoonDistance;
run PlotEarthMoon_system_frame.m
zvs=patch(isosurface(X,Y,Z,V2,0),'FaceColor','g','EdgeColor','none', 'FaceLighting', 'flat');
camlight; lighting phong
alpha(0.25)
view(46,28)
%--------------------------------------------------------------------------
%% Ploting zero velocity curves in 2D
[X, Y]=meshgrid(delta, delta);
rho1=sqrt((X+mu).^2 + Y.^2);
rho2=sqrt((X-(1-mu)).^2 + Y.^2);
V2=(X.^2+Y.^2)+2*(1-mu)./rho1+2*mu./rho2-C1;
UnitDistance=EarthMoonDistance;
run PlotEarthMoon_system_frame.m
contour(X,Y,V2,[0 0])