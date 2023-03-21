%%=========================================================================
%                        Cranfield University
%    Mathematics and Programming for Astrodynamics and Trajectory Design
%                      Student: Yi Qiang Ji Zhang
%                           
%        Copyright Cranfield University 2023, All Rights Reserved 
% 
%                
%%=========================================================================
%                            Assignment 1

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

% Spacecraft
spc_mass = 1250; % [kg] Dry mass
Isp = 310; % [s]
v_e = 9.81*Isp; % [m/s]

% Get astronomical parameters
[mu_Sun, mu_Earth, mu_Mars, R_Sun, R_Earth, R_Mars] = get_astro_constants();

% Planet number
n_Earth = 3; 
n_Mars = 4;

% Departure dates
date_dep_0_Earth = [2024, 01, 01, 0, 0, 0];
date_dep_1_Earth = [2025, 12, 31, 0, 0, 0];

% Convert to Julian Epoch
InTime_0 = date2mjd2000(date_dep_0_Earth);
InTime_1 = date2mjd2000(date_dep_1_Earth);

% Arrival dates
t_min = 200; % [days] Minimum transfer duration (aprox)
FinTime_0 = InTime_0 + t_min;
FinTime_1 = InTime_1 + t_min;

% Bounded values after anlaysis
InTime_0 = 8925;
InTime_1 = 9120;
FinTime_0 = 9225;
FinTime_1 = 9550;

% Departure and Arrival Grids
incr = 5;
DepartureGrid = InTime_0:incr:InTime_1;
ArrivalGrid = FinTime_0:incr:FinTime_1;

% Falcon 9
F9_C3_vec = [0 5 10 15 20 25 30 40];
F9_Wet_Mass_vec = [5750 5000 4250 3600 3000 2400 1900 1000];

% Initialize the delta_V solution matrix
DV_solutions = zeros(length(DepartureGrid), length(ArrivalGrid));
mass_solutions = zeros(length(DepartureGrid), length(DepartureGrid));
indexMin = zeros(length(DepartureGrid), length(ArrivalGrid));

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

% Max v_inf
v_inf_max = sqrt(F9_C3_vec(end)); % [km/s]


% Transfer delta_V using patched conics
for k=1:length(ArrivalGrid)
    for j=1:length(DepartureGrid)
    
        [r_dep_Earth,v_dep_Earth] = EphSS_car(n_Earth, DepartureGrid(j)); % [m and m/s]
        [r_arr_Mars,v_arr_Mars] = EphSS_car(n_Mars, ArrivalGrid(k)); % [m and m/s]
        
        delta_t_target = (ArrivalGrid(k) - DepartureGrid(j))*24*3600;
        
        % Short path
        [r1_dot_short, rf_dot_trans_real_short] = LambertArc(mu_Sun, r_dep_Earth, v_dep_Earth, ...
            r_arr_Mars, v_arr_Mars, t_m(1), delta_t_target);

        [delta_V_total_short] = get_total_deltaV(mu_Earth, mu_Mars, r_park, r_op, v_c_Earth, v_c_Mars, ...
            r1_dot_short, rf_dot_trans_real_short, v_dep_Earth, v_arr_Mars);

        % Long path
        [r1_dot_long, rf_dot_trans_real_long] = LambertArc(mu_Sun, r_dep_Earth, v_dep_Earth, ...
            r_arr_Mars, v_arr_Mars, t_m(2), delta_t_target);

        [delta_V_total_long] = get_total_deltaV(mu_Earth, mu_Mars, r_park, r_op, v_c_Earth, v_c_Mars, ...
            r1_dot_long, rf_dot_trans_real_long, v_dep_Earth, v_arr_Mars);

        % Delta_V matrix solutions
        [DV_solutions(j,k),indexMin(j,k)] = min([delta_V_total_short delta_V_total_long]);
        
        % Check which maneuvre is better
        if indexMin(j,k) == 1
            r1_dot = r1_dot_short;
            rf_dot = rf_dot_trans_real_short;
        else
            r1_dot = r1_dot_long;
            rf_dot = rf_dot_trans_real_long;
        end

        v_inf_Earth = norm(v_dep_Earth - r1_dot);
        v_inf_Mars = norm(v_arr_Mars - rf_dot);
        v_p_Mars = sqrt(mu_Mars*((2/r_op) + v_inf_Mars^2/mu_Mars));

        % Wet mass
        if v_inf_Earth <= v_inf_max
            C3 = v_inf_Earth^2;
            wet_mass = interp1(F9_C3_vec, F9_Wet_Mass_vec,C3);
        else
            v_inf_launcher = v_inf_max;
            C3 = v_inf_launcher^2;
            wet_mass = interp1(F9_C3_vec, F9_Wet_Mass_vec,C3, 'spline');
        end
        
        delta_V_2 = v_p_Mars - v_c_Mars;
        dry_mass = wet_mass*exp(-(delta_V_2*1000)/v_e);

        % Compute the dry mass matrix solutions
        mass_solutions(j,k) = dry_mass;
    end

end 

% Check boundaries of solutions
minDV = min(min(DV_solutions));
maxDV = max(max(DV_solutions));
minMass = min(min(mass_solutions));
maxMass = max(max(mass_solutions));

% Filter masses that are only >= spc_mass (dry mass)
mass_solutions_filtered = mass_solutions;
mass_solutions_filtered(mass_solutions_filtered <= spc_mass) = NaN;

%% Print results

% Find the indices of non NaN values of the matrix
[dep_dates, arr_dates] = find(~isnan(mass_solutions_filtered));

dep_dates_sol = zeros(length(dep_dates), 1);
arr_dates_sol = zeros(length(arr_dates), 1);

% Get departure and arrival dates
for i=1:length(dep_dates)
dep_dates_sol(i,1) = DepartureGrid(dep_dates(i));
arr_dates_sol(i,1) = ArrivalGrid(arr_dates(i));
end

% Find the delta_V and dry masses
delta_V_sol = DV_solutions(~isnan(DV_solutions));
mass_sol = mass_solutions_filtered(~isnan(mass_solutions_filtered));

dep_dates_sol_matrix = zeros(length(dep_dates), 6);

% Define the format string
formatSpec = '%-12s %-13s %-11s %-18s %-21s\n';
formatSpec2 = '%-12d %-13d %-11d %-18.2f %-21.2f\n';

% Define the column headers
headers = {'Dep year', 'Dep month', 'Dep day', 'Delta_V [km/s]', 'Dry mass [kg]'};
% Print the column headers
fprintf(formatSpec, headers{:});

% Print results
for i = 1:length(dep_dates)

    dep_dates_sol_aux = mjd20002date(dep_dates_sol(i,1));
    arr_dates_sol_aux = mjd20002date(arr_dates_sol(i,1));
    dep_dates_sol_matrix(i,:) = dep_dates_sol_aux;
    
    % Print results
    fprintf(formatSpec2, dep_dates_sol_matrix(i,1), ...
        dep_dates_sol_matrix(i,2), dep_dates_sol_matrix(i,3), delta_V_sol(i), mass_sol(i));
end

%% Plots

% Plot Pork-Chop plot
figure
hold on
contourf(DepartureGrid, ArrivalGrid, DV_solutions', minDV:0.1:8.5);
% plot(InTime_0, FinTime_1, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0, 0, 0], ...
%     'Marker', 'square', 'MarkerSize', 10, 'LineStyle', 'none');
col = colorbar;
col.Label.Interpreter = 'latex';
col.Label.String = '$\Delta V_{\mathrm{total}}$ [km/s]';
xlabel('Earth departure dates [MJD2000]')
ylabel('Mars arrival dates [MJD2000]')
title('\textbf{Pork Chop plot for Earth - Mars trajectory [Delta V]}')
hold off

figure
contourf(DepartureGrid, ArrivalGrid, mass_solutions',0:50:maxMass);
col = colorbar;
col.Label.Interpreter = 'latex';
col.Label.String = 'Mass [kg]';
xlabel('Earth departure dates [MJD2000]')
ylabel('Mars arrival dates [MJD2000]')
title('\textbf{Pork Chop plot for Earth - Mars trajectory [Mass]}')
hold off

figure
contourf(DepartureGrid, ArrivalGrid, mass_solutions',600:50:maxMass);
col = colorbar;
col.Label.Interpreter = 'latex';
col.Label.String = 'Mass [kg]';
xlabel('Earth departure dates [MJD2000]')
ylabel('Mars arrival dates [MJD2000]')
title('\textbf{Pork Chop plot for Earth - Mars trajectory [Mass semi-filtered]}')
hold off

figure
contourf(DepartureGrid, ArrivalGrid, mass_solutions_filtered',spc_mass:50:maxMass);
col = colorbar;
col.Label.Interpreter = 'latex';
col.Label.String = 'Mass [kg]';
xlabel('Earth departure dates [MJD2000]')
ylabel('Mars arrival dates [MJD2000]')
title('\textbf{Pork Chop plot for Earth - Mars trajectory [Valid Dry Masses]}')
hold off