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
%                         Exercise 4

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
Assessed exercise 1

%}


% Remark use only rf final and not all the iterations

%% Main

% Get astronomical parameters
[mu_Sun, mu_Earth, mu_Mars, R_Sun, R_Earth, R_Mars] = get_astro_constants();

% Planet number [1-9]
n_Earth = 3; 
n_Mars = 4;

% Departure window dates [Y M D h m s]
date_dep_0_Earth = [2024, 01, 01, 0, 0, 0];
date_dep_1_Earth = [2024, 01, 01, 0, 0, 0];


date_arr_Mars = [2025, 31, 12, 0, 0, 0];






% Time TGO Transfer
% delta_t_target = 215*24*3600; % [days]

% Convert to Julian Epoch
mjd2000_dep_Earth = date2mjd2000(date_dep_Earth);
mjd2000_arr_Mars = date2mjd2000(date_arr_Mars);

InTime = mjd2000_dep_Earth;
FinTime = mjd2000_arr_Mars;

% Generate Grid of conditions to analyse

% Range = 3*30; % [days]

DepartureGrid = InTime;
ArrivalGrid = FinTime;



% Departure

h0 = 250; % [km]
hf = 400; % [km]

% Circular parking orbit
r_park = R_Earth + h0; % [km]
v_c_Earth = sqrt(mu_Earth*(1/r_park)); % Velocity in parking orbit


% Mars Arrival
r_op = R_Mars + hf; % Operational orbit
v_c_Mars = sqrt(mu_Mars*(1/r_op)); % Velocity in operational orbit


% Transfer delta_V using patched conics

DV_solutions = zeros(length(DepartureGrid), length(ArrivalGrid));

% Short way in Lambert problem
t_m = [1 -1]; % Short and long path

    for k=1:length(ArrivalGrid)
    
        for j=1:length(DepartureGrid)
        
        [r_dep_Earth,v_dep_Earth] = EphSS_car(n_Earth, DepartureGrid(j)); % [m and m/s]
        [r_arr_Mars,v_arr_Mars] = EphSS_car(n_Mars, ArrivalGrid(k)); % [m and m/s]
        
        delta_t_target = (ArrivalGrid(k) - DepartureGrid(j))*24*3600;

        [r1_dot_short, rf_dot_trans_real_short] = LambertArc(mu, r_dep_Earth, v_dep_Earth, ...
            r_arr_Mars, v_arr_Mars, t_m(1), delta_t_target);

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

        disp

        % Save all the plots
%         DV_solutions(j,k) = min([delta_V_total_short]);

        end
        
    end

minDV = min(min(DV_solutions));


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

delta_t_min = sqrt(a_min^3/mu) * (pi - t_m*(beta_e - sin(beta_e))); % [s]

F = 1 - r2_norm/p_min*(1 - cos_theta);
G = (r2_norm*r1_norm)/(sqrt(mu*p_min)) * sin_theta;
        
r1_dot = G^(-1) * (r2 - F*r1);
delta_theta = acos(cos_theta);


end

% Algorith to compute Minimum Energy Transfer
function [rf] = FGKepler_trA(mu, r0, r0_dot, delta_theta)

r0_norm = norm(r0);

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

function [rf, rf_dot, delta_E_sol, F, G] = FGKepler_dt(mu, r0, r0_dot, delta_t)

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
rf_norm = norm(rf);

G_dot = 1 - a/rf_norm*(1 - cos(delta_E_sol));

rf_dot = 1/G*(G_dot*rf - r0);

end


function [STM] = STM_Lambert(mu, rf, delta_E, F, G, r0, r0_dot, delta_t)

r0_norm = norm(r0); 
r0_dot_norm = norm(r0_dot);
rf_norm = norm(rf);

a = mu/(2*mu/r0_norm - r0_dot_norm^2);

% Safety check

if a < 0
    warning ( 'sma seed error: Ihe shooting method did not converge ');
    r0_dot_norm = [NaN NaN NaN];

    return
end


F_dot = - (sqrt(mu*a))/(rf_norm*r0_norm) * sin(delta_E);
G_dot = 1 - a/rf_norm * (1 - cos(delta_E));

rf_dot = F_dot*r0 + G_dot*r0_dot;

C = a * sqrt(a^3/mu)*(3*sin(delta_E) - (2 + cos(delta_E))*delta_E) - ...
    a*delta_t*(1 - cos(delta_E));

delta_r = rf - r0;
delta_v = rf_dot - r0_dot;

STM = r0_norm/mu*(1 - F)*(delta_r'*r0_dot - delta_v'*r0) + ...
    C/mu*rf_dot'*r0_dot + G*eye(3);
end


function [r1_dot, rf_dot_trans_real] = LambertArc(mu, r_dep_Earth, v_dep_Earth, r_arr_Mars, v_arr_Mars, t_m, delta_t_target)

% Minimum Energy Transfer of the spacecraft
[~, ~, delta_t_min, r1_dot, ~] = MinETransfer(mu, r_dep_Earth, r_arr_Mars, t_m);

% Evaluate tolerance
tol = 1e-3; % [km]
error = 1e14; % Initialize error 

% Number of iterations
num_it = 10; % Number of iterations until reach the desired slution

% For loop of the Continuation Method
for it=1:num_it
    
    % Lambda
    Lambda = it/num_it;
    delta_t = Lambda*delta_t_target + (1 - Lambda)*delta_t_min;

    % Number of maximum iterations
    num_max_it = 25;
    num_iter = 0;
    error = 1e3; % Initiate error
    
    % Loop
    while (error >= tol) && (num_iter < num_max_it)
    a = mu/(2*mu/norm(r_dep_Earth) - norm(r1_dot)^2);

    if a < 0
        warning('sma seed error: The shooting method did not converge');
        r1_dot = [NaN NaN NaN];
        rf_dot_trans_real = [NaN NaN NaN];
        return 
    end

    [rf_trans_real, ~, delta_E, F, G] = FGKepler_dt(mu, r_dep_Earth, r1_dot, delta_t);
    error = norm(rf_trans_real - r_arr_Mars);
    
    [STM] = STM_Lambert(mu, rf_trans_real, delta_E, F, G, r_dep_Earth, r1_dot, delta_t);
    delta_r_dot_corr = inv(STM)*(rf_trans_real - r_arr_Mars)';
    
    r1_dot_new_trans = r1_dot - delta_r_dot_corr';
    r1_dot = r1_dot_new_trans;

    num_iter = num_iter + 1;    
    
    end

    if num_iter>num_max_it
        warning('num MaxIter issue : The shooting method did not converge')
        r1_dot = [NaN NaN NaN];
        rf_dot_trans_real = [NaN NaN NaN];
        return
    end

end

[rf_trans_real, rf_dot_trans_real, ~, ~, ~] = FGKepler_dt(mu, r_dep_Earth, r1_dot, delta_t_target);

end