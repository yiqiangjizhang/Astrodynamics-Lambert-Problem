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

% Problem statement
%{
Assessed Exercise 1. Your space agency is preparing a Mars remote sensing observation
mission. The mission is intended to be launched between 1/01/2024 and 31/12/2025. The
orbiter has a mass of 1250 kg and a propulsion system with an Isp of 310s. The targeted
final orbit is a 400 km circular orbit.

Perform a complete assessment of Earth-Mars direct insertion launch opportunities for a
Falcon 9 Heavy launch with booster recovery (see below performance of this launcher). Your
space agency programme manager is interested to know the extend of the launch window
opportunity and your recommendation as baseline launch date.

Note: Since the objective is to insert the spacecraft in its operational circular orbit, a launch
opportunity will be defined as any launch date that allows you to insert the spacecraft dry
mass into the 400 km circular orbit. The inclination of the Mars circular orbit is not
constrained (which is to say that is ignored at this stage).

Hence, the objective is to construct a script or a function that given a departure date and a
time of flight (or an arrival date) it computes how much mass Falcon 9 inserts into
interplanetary orbit. Next, with this mass and the Isp, the final mass into the final Mars orbit
can be computed.

Note that Falcon 9 launch performance is provided in slide 29. Neglect launch adapter mass.
%}

%% Setup workspace

% Clear workspace, command window and close windows
clear;
close all;
clc;

% Set interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% Main

% Spacecraft
spc_mass = 1250; % [kg] Dry mass
Isp = 310; % [s]
v_e = 9.81 * Isp; % [m/s]

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
t_min = 215; % [days] Minimum transfer duration (aprox)
FinTime_0 = InTime_0 + t_min;
FinTime_1 = InTime_1 + t_min;

% Bounded values after analysis
InTime_0 = 8925;
InTime_1 = 9120;
FinTime_0 = 9225;
FinTime_1 = 9550;

% Departure and Arrival Grids
incr = 1; % [days]
DepartureGrid = InTime_0:incr:InTime_1;
ArrivalGrid = FinTime_0:incr:FinTime_1;

% Falcon 9 specs
F9_C3_vec = [0 5 10 15 20 25 30 40];
F9_Wet_Mass_vec = [5750 5000 4250 3600 3000 2400 1900 1000];

% Initialize the delta_V solution matrix
DV_solutions = zeros(length(DepartureGrid), length(ArrivalGrid));
mass_solutions = zeros(length(DepartureGrid), length(DepartureGrid));
indexMin = zeros(length(DepartureGrid), length(ArrivalGrid));

% Plot matrices
C3_matrix = zeros(length(DepartureGrid), length(DepartureGrid));
TOF_matrix = zeros(length(DepartureGrid), length(DepartureGrid));
v_inf_matrix = zeros(length(DepartureGrid), length(DepartureGrid));

% Transfer method (short or long path)
t_m = [1 -1]; % Short and long path respectively

% Departure and Arrival altitudes
h0 = 250; % [km]
hf = 400; % [km]

% Earth circular parking orbit
r_park = R_Earth + h0; % [km]
v_c_Earth = sqrt(mu_Earth * (1 / r_park)); % Velocity in parking orbit

% Mars Arrival
r_op = R_Mars + hf; % Operational orbit
v_c_Mars = sqrt(mu_Mars * (1 / r_op)); % Velocity in operational orbit

% Max v_inf
v_inf_max = sqrt(F9_C3_vec(end)); % [km/s]

% Transfer delta_V using patched conics
for k = 1:length(ArrivalGrid)

    for j = 1:length(DepartureGrid)

        % Ephemerides of Earth and Mars
        [r_dep_Earth, v_dep_Earth] = EphSS_car(n_Earth, DepartureGrid(j)); % [m and m/s]
        [r_arr_Mars, v_arr_Mars] = EphSS_car(n_Mars, ArrivalGrid(k)); % [m and m/s]

        % Time of Flight
        delta_t_target = (ArrivalGrid(k) - DepartureGrid(j)) * 24 * 3600; % [s]

        % Short path
        [r1_dot_short, rf_dot_trans_real_short] = LambertArc(mu_Sun, r_dep_Earth, ...
            r_arr_Mars, t_m(1), delta_t_target);

        [delta_V_total_short] = get_total_deltaV(mu_Earth, mu_Mars, r_park, r_op, v_c_Earth, v_c_Mars, ...
            r1_dot_short, rf_dot_trans_real_short, v_dep_Earth, v_arr_Mars);

        % Long path
        [r1_dot_long, rf_dot_trans_real_long] = LambertArc(mu_Sun, r_dep_Earth, ...
            r_arr_Mars, t_m(2), delta_t_target);

        [delta_V_total_long] = get_total_deltaV(mu_Earth, mu_Mars, r_park, r_op, v_c_Earth, v_c_Mars, ...
            r1_dot_long, rf_dot_trans_real_long, v_dep_Earth, v_arr_Mars);

        % Delta_V matrix solutions
        [DV_solutions(j, k), indexMin(j, k)] = min([delta_V_total_short delta_V_total_long]);

        % Check which maneuvre is better
        if indexMin(j, k) == 1
            r1_dot = r1_dot_short;
            rf_dot = rf_dot_trans_real_short;
        else
            r1_dot = r1_dot_long;
            rf_dot = rf_dot_trans_real_long;
        end

        % Compute hyperbolic v_inf at Earth and Mars
        v_inf_Earth = norm(v_dep_Earth - r1_dot);
        v_inf_Mars = norm(v_arr_Mars - rf_dot);
        v_p_Mars = sqrt(mu_Mars * ((2 / r_op) + v_inf_Mars ^ 2 / mu_Mars));

        % Compute wet mass
        if v_inf_Earth <= v_inf_max
            C3_matrix(j, k) = v_inf_Earth ^ 2;
            wet_mass = interp1(F9_C3_vec, F9_Wet_Mass_vec, C3_matrix(j, k));
        else
            v_inf_launcher = v_inf_max;
            C3_matrix(j, k) = v_inf_launcher ^ 2;
            wet_mass = interp1(F9_C3_vec, F9_Wet_Mass_vec, C3_matrix(j, k), 'spline');
        end

        % Compute required delta_V at arrival to Mars
        delta_V_2 = v_p_Mars - v_c_Mars;

        % Compute the dry mass matrix solutions
        dry_mass = wet_mass * exp(- (delta_V_2 * 1000) / v_e);
        mass_solutions(j, k) = dry_mass;

        % TOF days matrix
        TOF_matrix(j, k) = ArrivalGrid(k) - DepartureGrid(j); % [days]

        % Arrival v_inf matrix
        v_inf_matrix(j, k) = v_inf_Mars; % [km/s]
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
for i = 1:length(dep_dates)
    dep_dates_sol(i, 1) = DepartureGrid(dep_dates(i));
    arr_dates_sol(i, 1) = ArrivalGrid(arr_dates(i));
end

% Find the delta_V and dry masses
delta_V_sol = DV_solutions(~isnan(DV_solutions));
mass_sol = mass_solutions_filtered(~isnan(mass_solutions_filtered));

dep_dates_sol_matrix = zeros(length(dep_dates), 6);

% Final solution matrix
launch_parameters = zeros(length(dep_dates), 5);

% Define the format string
formatSpec = '%-12s %-13s %-11s %-18s %-21s\n';
formatSpec2 = '%-12d %-13d %-11d %-18.2f %-21.2f\n';

% Define the column headers
headers = {'Dep year', 'Dep month', 'Dep day', 'Delta_V [km/s]', 'Dry mass [kg]'};
% Print the column headers
fprintf(formatSpec, headers{:});

% Print results
for i = 1:length(dep_dates)

    % Convert MJD2000 to date
    dep_dates_sol_aux = mjd20002date(dep_dates_sol(i, 1));
    arr_dates_sol_aux = mjd20002date(arr_dates_sol(i, 1));
    dep_dates_sol_matrix(i, :) = dep_dates_sol_aux;

    % Fill the final solution matrix
    launch_parameters(i, 1) = dep_dates_sol_matrix(i, 1);
    launch_parameters(i, 2) = dep_dates_sol_matrix(i, 2);
    launch_parameters(i, 3) = dep_dates_sol_matrix(i, 3);
    launch_parameters(i, 4) = delta_V_sol(i);
    launch_parameters(i, 5) = mass_sol(i);

    % Print results
    fprintf(formatSpec2, dep_dates_sol_matrix(i, 1), ...
        dep_dates_sol_matrix(i, 2), dep_dates_sol_matrix(i, 3), delta_V_sol(i), mass_sol(i));
end

%% Plots

% Plot Delta V Pork-Chop plot
figure('position', [500, 500, 1000, 420])
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

% Plot Dry mass Pork-Chop plot
figure('position', [500, 500, 1000, 420])
contourf(DepartureGrid, ArrivalGrid, mass_solutions', 0:50:maxMass);
col = colorbar;
col.Label.Interpreter = 'latex';
col.Label.String = 'Mass [kg]';
xlabel('Earth departure dates [MJD2000]')
ylabel('Mars arrival dates [MJD2000]')
title('\textbf{Pork Chop plot for Earth - Mars trajectory [Mass]}')
hold off

% Plot Dry mass Pork-Chop plot from 600 kg to max
figure('position', [500, 500, 1000, 420])
contourf(DepartureGrid, ArrivalGrid, mass_solutions', 600:50:maxMass);
col = colorbar;
col.Label.Interpreter = 'latex';
col.Label.String = 'Mass [kg]';
xlabel('Earth departure dates [MJD2000]')
ylabel('Mars arrival dates [MJD2000]')
title('\textbf{Pork Chop plot for Earth - Mars trajectory [Mass semi-filtered]}')
hold off

% Plot Dry mass Pork-Chop plot that are valid (>=1250 kg, spacecraft dry mass)
figure('position', [500, 500, 1000, 420])
contourf(DepartureGrid, ArrivalGrid, mass_solutions_filtered', spc_mass:50:maxMass);
col = colorbar;
col.Label.Interpreter = 'latex';
col.Label.String = 'Mass [kg]';
xlabel('Earth departure dates [MJD2000]')
ylabel('Mars arrival dates [MJD2000]')
title('\textbf{Pork Chop plot for Earth - Mars trajectory [Valid Dry Masses]}')
hold off

% Overall pork chop plot with launcher C3 parameter
figure('position', [500, 500, 1000, 420])

% Set levels to plot
C3_levels = 11:2:40;
TOF_levels = 100:50:600;
v_inf_levels = 1:1:10;

% Colors
col1 = [0.8, 0.2, 0.2];
col2 = [0.2, 0.2, 0.8];
col3 = [0.4, 0.4, 0.4];
col4 = [0.9290 0.6940 0.1250];

% Plot
set(gcf, 'color', 'w')
hold on
box on
grid on
grid minor
[c1, h1] = contour(DepartureGrid, ArrivalGrid, v_inf_matrix', v_inf_levels, 'color', col1, 'linewidth', 1.25);
[c2, h2] = contour(DepartureGrid, ArrivalGrid, C3_matrix', C3_levels, 'color', col2, 'linewidth', 1.25);
[c3, h3] = contour(DepartureGrid, ArrivalGrid, TOF_matrix', TOF_levels, 'color', col3);
[c4, h4] = contour(DepartureGrid, ArrivalGrid, DV_solutions', 5.5:0.5:8.5, 'color', col4, 'linewidth', 1.25);
clabel(c1, h1, 'Color', col1)
clabel(c2, h2, 'Color', col2)
clabel(c3, h3, 'Color', col3)
clabel(c4, h4, 'Color', col3)
xlabel('Earth departure dates [MJD2000]')
ylabel('Mars arrival dates [MJD2000]')
title('\textbf{Pork Chop plot for Earth - Mars trajectory [Valid Dry Masses]}')
legend({'$v_{\infty \, \mathrm{arrival}}$ [km/s]', '${\rm C}_3$ Proton [km$^2$/s$^2$]', 'TOF [days]', '$\Delta V_{\mathrm{total}}$ [km/s]'}, 'Location', 'northeastoutside')
hold off
