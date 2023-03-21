%%=========================================================================
%                        Cranfield University
%    Mathematics and Programming for Astrodynamics and Trajectory Design
%                  Student: Yi Qiang Ji Zhang
%                           
%     Copyright Cranfield University 2021, All Rights Reserved 
% 
% Versions:
%           v0.1 file created
%                
%%=========================================================================
%                         Exercise 1

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
Exercise 1. ExoMars Trace Gas Orbiter departed Earth on 14/03/2016, and
arrived at Mars on 15/10/2016. What is the minimum energy orbit that would 
link the position of Earth and Mars those two days? What is the time of 
flight of the minimum energy orbit? Was this the trajectory followed 
by ExoMars TGO? 
%}


%% Main Exercise 1

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
[a_min, e_min, delta_t_min] = MinETransfer(r_dep_Earth, r_arr_Mars, t_m);

disp("Minimum time: " + delta_t_min + " days")


%% Functions

% Algorith to compute Minimum Energy Transfer
function [a_min,e_min, delta_t_min] = MinETransfer(r1, r2, t_m)

% MinETransfer.m - Calculate the Minimum Energy Transfer (Lambert Arch).
%
% PROTOTYPE:
%   [a_min, e_min, delta_t_min] = MinETransfer(r_dep_Earth, r_arr_Mars, t_m);
%
% DESCRIPTION:
%   Computes the Minimum Energy Transfer using general approach to
%   Lambert's problem
%
% INPUT:
%   a_min[1] 
%   r[1,3]   Cartesian position of the body (Sun-centered for all bodies,
%            Earth-centered for the Moon) [km].
%   v[1,3]   Cartesian velocity [km/s].
%   t_m[1]   Transfer Method. Value is either "+1" for Short Way or "-1" for
%            Long Way
%
% OUTPUT:
%   a_min[1] Minimum semi-major axis [km]
%   e_min[1] Minimum eccentricity [km]
%   delta_t_min[1] Minimum Transfer time [days]
%   r[1,3]   Cartesian position of the body (Sun-centered for all bodies,
%            Earth-centered for the Moon).
%   v[1,3]   Cartesian velocity.
%

% Sun parameter
mu = getAstroConstants('Sun', 'mu'); % [km3/s2]

c = norm(r2 - r1);
r1_norm = norm(r1);
r2_norm = norm(r2);

a_min = 0.25*(r1_norm + r2_norm + c);
cos_theta = (dot(r1,r2))/(r1_norm*r2_norm);
p_min = (r1_norm*r2_norm)/c * (1 - cos_theta);
e_min = sqrt(1 - p_min/a_min);
beta_e = 2 * asin(sqrt((2*a_min - c)/(2*a_min)));

delta_t_min = sqrt(a_min^3/mu) * (pi - t_m*(beta_e - sin(beta_e)));

% Convert to days
delta_t_min = delta_t_min/(24*3600);

end


