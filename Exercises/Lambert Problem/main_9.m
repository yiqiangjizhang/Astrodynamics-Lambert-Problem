%%=========================================================================
%                        Cranfield University
%    Mathematics and Programming for Astrodynamics and Trajectory Design
%                      Student: Yi Qiang Ji Zhang
%                           
%        Copyright Cranfield University 2023, All Rights Reserved 
% 
%                
%%=========================================================================
%                            Exercise 6

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
Exercise 9. 
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
InTime = date2mjd2000(date_dep_Earth);
FinTime = date2mjd2000(date_arr_Mars);

% Generate Grid of conditions to analyse
Range = 1; % [days]
DepartureGrid = InTime - Range:5:InTime + Range;
ArrivalGrid = FinTime - Range:5:FinTime + Range;

% Initialize the delta_V solution matrix
DV_solutions = zeros(length(DepartureGrid), length(ArrivalGrid));

% Short way in Lambert problem
t_m = [1 -1]; % Short and long path

% Departure and Arrival altitudes
h0 = 250; % [km]
hf = 400; % [km]

% Earth circular parking orbit
r_park = R_Earth + h0; % [km]
v_c_Earth = sqrt(mu_Earth*(1/r_park)); % Velocity in parking orbit

% Mars Arrival
r_op = R_Mars + hf; % Operational orbit
v_c_Mars = sqrt(mu_Mars*(1/r_op)); % Velocity in operational orbit


% Transfer delta_V using patched conics
for k=1:length(ArrivalGrid)
    for j=1:length(DepartureGrid)
    
        [r_dep_Earth,v_dep_Earth] = EphSS_car(n_Earth, DepartureGrid(j)); % [m and m/s]
        [r_arr_Mars,v_arr_Mars] = EphSS_car(n_Mars, ArrivalGrid(k)); % [m and m/s]
        
        delta_t_target = (ArrivalGrid(k) - DepartureGrid(j))*24*3600;
        
        [r1_dot_short, rf_dot_trans_real_short] = LambertArc(mu_Sun, r_dep_Earth, ...
            r_arr_Mars, t_m(1), delta_t_target);
        
        [r1_dot_long, rf_dot_trans_real_long] = LambertArc(mu_Sun, r_dep_Earth, ...
            r_arr_Mars, t_m(2), delta_t_target);
        
        % Delta V Short
        v_inf_dep = norm(r1_dot_short - v_dep_Earth);
        v_inf_arr = norm(rf_dot_trans_real_short - v_arr_Mars);

        % Departure
        v_p_Earth = sqrt(mu_Earth*(2/r_park + v_inf_dep^2/mu_Earth)); % Velocity at periapsis in scape orbit
        delta_v_1 = v_p_Earth - v_c_Earth;
        
        % Arrival
        v_p_Mars = sqrt(mu_Mars*(2/r_op + v_inf_arr^2/mu_Mars)); % Velocity at periapsis in scape orbit
        delta_v_2 = v_p_Mars - v_c_Mars;

        % Total delta V
        delta_V_total_short = delta_v_1 + delta_v_2;
        
        % Display total delta V
        disp("Total delta V = " + delta_V_total_short + " [km/s]");
        
    end

end

