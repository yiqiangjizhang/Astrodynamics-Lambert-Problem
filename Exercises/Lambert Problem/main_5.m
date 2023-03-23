%%=========================================================================
%                        Cranfield University
%    Mathematics and Programming for Astrodynamics and Trajectory Design
%                      Student: Yi Qiang Ji Zhang
%                           
%        Copyright Cranfield University 2023, All Rights Reserved 
% 
%                
%%=========================================================================
%                            Exercise 5

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
Exercise 5. 
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

% Time TGO Transfer
delta_t = 215*24*3600; % [days]

% Tolerance and error
tol = 1e-3; % [km]
error = 1e9; 

% Number of iterations
it = 1;

% Minimize error
while error >= tol

    % Final position
    [rf_trans_real, ~, delta_E_sol, F, G] = FGKepler_dt(mu_Sun, r_dep_Earth, r1_dot, delta_t);
    
    % Compute error
    error = norm(rf_trans_real - r_arr_Mars);
    
    % Correct error using State Transition Matrix
    [STM] = STM_Lambert(mu_Sun, rf_trans_real, delta_E_sol, F, G, r_dep_Earth, r1_dot, delta_t);
    delta_r_dot_corr = STM\(rf_trans_real - r_arr_Mars)';
    r1_dot_new_trans = r1_dot - delta_r_dot_corr';
    
    % Compute corrected initial velocity
    r1_dot = r1_dot_new_trans;
    
    it = it + 1;
end

% Error
error = norm(rf_trans_real - r_arr_Mars);

% Display answer
disp("Error: " + error);