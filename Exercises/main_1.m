%%=========================================================================
%                        Cranfield University
%    Mathematics and Programming for Astrodynamics and Trajectory Design
%                      Student: Yi Qiang Ji Zhang
%                           
%        Copyright Cranfield University 2023, All Rights Reserved 
% 
%                
%%=========================================================================
%                            Exercise 1

%% Setup workspace

% Clear workspace, command window and close windows
clear;
close all;
clc;

% Set interpreter to latex
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


%% Problem statement

%{
Exercise 1. ExoMars Trace Gas Orbiter departed Earth on 14/03/2016, and
arrived at Mars on 15/10/2016. What is the minimum energy orbit that would 
link the position of Earth and Mars those two days? What is the time of 
flight of the minimum energy orbit? Was this the trajectory followed 
by ExoMars TGO? 
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
t_m = -1;

% Minimum Energy Transfer
[a_min, e_min, delta_t_min] = MinETransfer(mu_Sun, r_dep_Earth, r_arr_Mars, t_m);

% Convert to days
delta_t_min = delta_t_min/(24*3600);

% Print Minimum time [days]
disp("Minimum time: " + delta_t_min + " days");