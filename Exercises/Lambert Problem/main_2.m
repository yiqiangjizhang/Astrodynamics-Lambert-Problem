%%=========================================================================
%                        Cranfield University
%    Mathematics and Programming for Astrodynamics and Trajectory Design
%                      Student: Yi Qiang Ji Zhang
%                           
%        Copyright Cranfield University 2023, All Rights Reserved 
% 
%                
%%=========================================================================
%                            Exercise 2

%% Setup workspace

% Clear workspace, command window and close windows
clear;
close all;
clc;

% Set interpreter to latex
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


% Problem statement
%{
Exercise 2. ExoMars Trace Gas Orbiter departed Earth on 14/03/2016, and arrived at
Mars on 15/10/2016.
Program the F and G solutions to the two body problem. Verify the answer
by comparing it to a numerical integration of the differential equations of
motion:

r''+ mu_E/r3 r = 0

Using the F and G techniques, plot the orbit of Mars, Earth and the
minimum energy transfer for ExoMars TGO.

%}

%% Main

% Get astronomical parameters
[mu_Sun, mu_Earth, mu_Mars, R_Sun, R_Earth, R_Mars] = get_astro_constants();

% Planet number
n_Earth = 3; 
n_Mars = 4;

% Departure and arrival dates
date_dep_Earth = [2016, 3, 14, 0, 0, 0];
date_arr_Mars = [2016, 10, 15, 0, 0, 0];

% Convert to Julian Epoch
mjd2000_dep_Earth = date2mjd2000(date_dep_Earth);
mjd2000_arr_Mars = date2mjd2000(date_arr_Mars);

% Get ephemeris of Earth and Mars
[r_dep_Earth,v_dep_Earth] = EphSS_car(n_Earth,mjd2000_dep_Earth); % [m and m/s]
[r_arr_Mars,v_arr_Mars] = EphSS_car(n_Mars,mjd2000_arr_Mars); % [m and m/s]

% Short path in Lambert problem
t_m = 1;

% Minimum Energy Transfer
[a_min, e_min, delta_t_min, r1_dot, delta_theta] = MinETransfer(mu_Sun, r_dep_Earth, r_arr_Mars, t_m);

% Propagate orbits

% Full propagation
full_delta_theta = 0:0.01:2*pi; % Planet's propagation
propagation_theta = linspace(0,delta_theta,1000); % Transfer propagation

% Propagation of Earth, Transfer Orbit and Mars Orbit
[rf_Earth] = FGKepler_trA(mu_Sun, r_dep_Earth, v_dep_Earth, full_delta_theta);
[rf_trans] = FGKepler_trA(mu_Sun, r_dep_Earth, r1_dot, propagation_theta);
[rf_Mars] = FGKepler_trA(mu_Sun, r_arr_Mars, v_arr_Mars, full_delta_theta);

% Answer
disp("r1_dot: ")
disp(r1_dot)
disp("[km/s]")

% Plot propagation
figure
plot3(rf_Earth(1,:), rf_Earth(2,:), rf_Earth(3,:));
hold on
plot3(rf_trans(1,:), rf_trans(2,:), rf_trans(3,:));
plot3(rf_Mars(1,:), rf_Mars(2,:), rf_Mars(3,:));
xlabel('x')
ylabel('y')
title('\textbf{Transfer Orbit}')
grid on
grid minor
legend('Earth Orbit', 'Transfer Orbit', 'Mars Orbit', 'Sun')
HFIG=drawPlanet('Sun',[0 0 0],1,25);
% HFIG2=drawPlanet('Earth',r_dep_Earth,1,100);
% HFIG3=drawPlanet('Mars',r_arr_Mars,1,100);
hold off
