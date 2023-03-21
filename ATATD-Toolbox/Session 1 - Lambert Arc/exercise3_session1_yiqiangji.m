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
%                         Exercise 3

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
Exercise 3. ExoMars Trace Gas Orbiter departed Earth on 14/03/2016, and arrived at
Mars on 15/10/2016.
Program the F and G solutions to the two body problem. Verify the answer
by comparing it to a numerical integration of the differential equations of
motion:

r''+ mu_E/r3 r = 0

Using the F and G techniques, plot the orbit of Mars, Earth and the
minimum energy transfer for ExoMars TGO.

%}


% Remark use only rf final and not all the iterations

%% Main


% Sun parameter
mu = getAstroConstants('Sun', 'mu'); % [km3/s2]

% Get Ephemeris for Earth

% Planet number
n_Earth = 3; 
n_Mars = 4;

% Departure and arrival dates
date_dep_Earth = [2016, 3, 14, 0, 0, 0];
date_arr_Mars = [2016, 10, 15, 0, 0, 0];

% Convert to Julian Epoch
mjd2000_dep_Earth = date2mjd2000(date_dep_Earth);
mjd2000_arr_Mars = date2mjd2000(date_arr_Mars);

[r_dep_Earth,v_dep_Earth] = EphSS_car(n_Earth,mjd2000_dep_Earth); % [m and m/s]
[r_arr_Mars,v_arr_Mars] = EphSS_car(n_Mars,mjd2000_arr_Mars); % [m and m/s]

% Short way in Lambert problem
t_m = 1;

% Minimum Energy Transfer of the spacecraft
[a_min, e_min, delta_t_min, r1_dot, delta_theta] = MinETransfer(mu, r_dep_Earth, r_arr_Mars, t_m);

% Propagate orbits

% Full propagation
full_delta_theta = 0:0.01:2*pi; % Planet's propagation
propagation_theta = linspace(0,delta_theta,1000); % Transfer propagation

% Propagation of Earth, Transfer Orbit and Mars Orbit
[rf_Earth] = FGKepler_trA(r_dep_Earth, v_dep_Earth, full_delta_theta);
[rf_trans] = FGKepler_trA(r_dep_Earth, r1_dot, propagation_theta);
[rf_Mars] = FGKepler_trA(r_arr_Mars, v_arr_Mars, full_delta_theta);

% Time TGO Transfer
delta_t = 215*24*3600; % [days]
[rf_trans_real] = FGKepler_dt(r_dep_Earth, r1_dot, delta_t);


% Should use
% error = norm(rf_trans_real - rf_trans(:,end));

% No problem, it's pathed conics
error = norm(rf_trans_real - r_arr_Mars);

disp("Error: " + error);


%% Functions

% Algorith to compute Minimum Energy Transfer
function [a_min, e_min, delta_t_min, r1_dot, delta_theta] = MinETransfer(mu, r1, r2, t_m)

c = norm(r2 - r1);
r1_norm = norm(r1);
r2_norm = norm(r2);

a_min = 0.25*(r1_norm + r2_norm + c);

cos_theta = (dot(r1,r2))/(r1_norm*r2_norm);
sin_theta = t_m * sqrt(1 - cos_theta^2);

p_min = (r1_norm*r2_norm)/c * (1 - cos_theta);
e_min = sqrt(1 - p_min/a_min);
beta_e = 2 * asin(sqrt((2*a_min - c)/(2*a_min)));

delta_t_min = sqrt(a_min^3/mu) * (pi - t_m*(beta_e - sin(beta_e)));

% Convert to days
delta_t_min = delta_t_min/(24*3600);

F = 1 - r2_norm/p_min*(1 - cos_theta);
G = (r2_norm*r1_norm)/(sqrt(mu*p_min)) * sin_theta;
        
r1_dot = G^(-1) * (r2 - F*r1);
delta_theta = acos(cos_theta);


end

% Algorith to compute Minimum Energy Transfer
function [rf] = FGKepler_trA(r0, r0_dot, delta_theta)

r0_norm = norm(r0);
mu = getAstroConstants('Sun', 'mu'); % [km3/s2]

% Angular momentum
h_vec = cross(r0, r0_dot);
h = norm(h_vec);
p = h^2/mu;

sigma0 = dot(r0,r0_dot)/sqrt(mu);

rf = zeros(3,length(delta_theta));

for i = 1:length(delta_theta)
    rf_norm = p*r0_norm/(r0_norm + (p - r0_norm)*cos(delta_theta(i)) - sqrt(p) * sigma0 * sin(delta_theta(i)));
    
    F = 1 - rf_norm/p * (1 - cos(delta_theta(i)));
    G = (rf_norm*r0_norm)/sqrt(mu*p) * sin(delta_theta(i));
    
    rf_aux = F*r0 + G*r0_dot;
    rf(1,i) = rf_aux(1);
    rf(2,i) = rf_aux(2);
    rf(3,i) = rf_aux(3);
end

end



function [rf] = FGKepler_dt(r0, r0_dot, delta_t)

% Sun parameter
mu = getAstroConstants('Sun', 'mu'); % [km3/s2]

r0_norm = norm(r0);
r0_dot_norm = norm(r0_dot);

a = mu/(2*mu/r0_norm - r0_dot_norm^2);

n = sqrt(mu/a^3);

delta_M = n*delta_t;

sigma_0 = (dot(r0,r0_dot))/sqrt(mu);

fun = @(delta_E) delta_E - (1 - r0_norm/a)*sin(delta_E) - sigma_0/sqrt(a) * (cos(delta_E)- 1) - delta_M; % function
x0 = delta_M; % initial point
delta_E_sol = fzero(fun,x0);

% Check
fun(delta_E_sol);

F = 1 - a/r0_norm * (1 - cos(delta_E_sol));
G = delta_t + sqrt(a^3/mu) * (sin(delta_E_sol) - delta_E_sol);
rf = F*r0 + G*r0_dot;


end
