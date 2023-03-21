%%=========================================================================
%                        Cranfield University
%    Mathematics and Programming for Astrodynamics and Trajectory Design
%                  Student: Yi Qiang Ji Zhang
%                           
%       Copyright Cranfield University 2021, All Rights Reserved 
% 
%%=========================================================================
%                         Exercise 0
%                    Name

clear all
close all
clc


%% Main

% Data
M_Earth = 5.972E24; % [kg]
G = 6.6743E-11; % [m3/kg/s2]
mu = G*M_Earth; % [m3/s2]


% Initial conditions
SV0 = [-18.676 6,246 12,474 0.551 -1.946 -3.886]; % [km km km km/s km/s km/s]

% Time span
tspan = 12*3600; % [s]

% ODE Equation r'' + mu/r3 *
F2BDyn = @(t,x) [x(4); % [dVx/dt]
    x(5); % [dVx/dt]
    x(6); % [dVx/dt]
    -muEarth*x(1)/(sqrt(x(1)^2 + x(2)^2 + x(3)^2)^3); % [dVx/dt]
    -muEarth*x(2)/(sqrt(x(1)^2 + x(2)^2 + x(3)^2)^3); % [dVx/dt]
    -muEarth*x(3)/(sqrt(x(1)^2 + x(2)^2 + x(3)^2)^3)]; % [dVx/dt]

% Propagations
options=odeset ('RelTol', 1e-6, 'AbsTol' , 1e-6);
[Tstep, SVt]=ode45(F2BDyn, [0, tf], SV0, options);


figure
plot3 (SVt (:,1) , SVt (:,2) , SVt(:,3));


Re = 6378; % [km]
figure
axis equal
hold on
[X, Y, Z] = sphere(20);
X = X*Re;
Y = Y*Re;
Z = Z*Re;
surf (X, Y, Z)
plot3 (SVt(:,1), SVt(:,2), SVt(:,3));
view (75, 20)

