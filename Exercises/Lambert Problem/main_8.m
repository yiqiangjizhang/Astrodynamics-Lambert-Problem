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
Exercise 3. 
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
Range = 3*30; % [days]
DepartureGrid = InTime - Range:5:InTime + Range;
ArrivalGrid = FinTime - Range:5:FinTime + Range;

% Initialize the delta_V solution matrix
DV_solutions = zeros(length(DepartureGrid), length(ArrivalGrid));

% Short way in Lambert problem
t_m = [1 -1]; % Short and long path


for k=1:length(ArrivalGrid)
    for j=1:length(DepartureGrid)
    
        [r_dep_Earth,v_dep_Earth] = EphSS_car(n_Earth, DepartureGrid(j)); % [m and m/s]
        [r_arr_Mars,v_arr_Mars] = EphSS_car(n_Mars, ArrivalGrid(k)); % [m and m/s]
        
        delta_t_target = (ArrivalGrid(k) - DepartureGrid(j))*24*3600;
        
        [r1_dot_short, rf_dot_trans_real_short] = LambertArc(mu_Sun, r_dep_Earth, v_dep_Earth, ...
            r_arr_Mars, v_arr_Mars, t_m(1), delta_t_target);
        
        [r1_dot_long, rf_dot_trans_real_long] = LambertArc(mu_Sun, r_dep_Earth, v_dep_Earth, ...
            r_arr_Mars, v_arr_Mars, t_m(2), delta_t_target);
        
        % Delta V Short
        delta_V_dep_trans_Earth_short = norm(r1_dot_short - v_dep_Earth);
        delta_V_arr_trans_Mars_short = norm(rf_dot_trans_real_short - v_arr_Mars);
        
        % Total delta V
        delta_V_total_short = delta_V_dep_trans_Earth_short + delta_V_arr_trans_Mars_short;
        
        % Delta V Long
        delta_V_dep_trans_Earth_long = norm(r1_dot_long - v_dep_Earth);
        delta_V_arr_trans_Mars_long = norm(rf_dot_trans_real_long - v_arr_Mars);
        
        % Total delta V Long
        delta_V_total_long = delta_V_dep_trans_Earth_long + delta_V_arr_trans_Mars_long;
      
        % Save all the plots
        DV_solutions(j,k) = min([delta_V_total_short delta_V_total_long]);
    end

end

minDV = min(min(DV_solutions));

% Plot Pork-Chop plot
figure
hold on
contourf(DepartureGrid, ArrivalGrid, DV_solutions', minDV:0.1:8.5);
plot(InTime, FinTime, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0, 0, 0], ...
    'Marker', 'square', 'MarkerSize', 10, 'LineStyle', 'none');
col = colorbar;
col.Label.String = '$\Delta V_{\mathrm{total}}$ [km/s]';
col.Label.Interpreter = 'latex';
xlabel('Earth departure dates [MJD2000]')
ylabel('Mars arrival dates [MJD2000]')
title('\textbf{Pork Chop plot for Earth - Mars trajectory }')