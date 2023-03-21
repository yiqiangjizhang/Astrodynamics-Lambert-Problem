
function [r1_dot, rf_dot_trans_real] = LambertArc(mu, r_dep_Earth, v_dep_Earth, r_arr_Mars, v_arr_Mars, t_m, delta_t_target)

% Minimum Energy Transfer of the spacecraft
[~, ~, delta_t_min, r1_dot, ~] = MinETransfer(mu, r_dep_Earth, r_arr_Mars, t_m);

% Evaluate tolerance
tol = 1e-3; % [km]

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
    error = 1e14; % Initialize error 
    
    % Loop
    while (error >= tol) && (num_iter < num_max_it)
    a = mu/(2*mu/norm(r_dep_Earth) - norm(r1_dot)^2);

    % Check sma (semi-major axis is negative)
    if a < 0
        warning('a seed error: The shooting method did not converge');
        r1_dot = [NaN NaN NaN];
        rf_dot_trans_real = [NaN NaN NaN];
        return 
    end

    [rf_trans_real, ~, delta_E, F, G] = FGKepler_dt(mu, r_dep_Earth, r1_dot, delta_t);
    error = norm(rf_trans_real - r_arr_Mars);
    
    [STM] = STM_Lambert(mu, rf_trans_real, delta_E, F, G, r_dep_Earth, r1_dot, delta_t);
    delta_r_dot_corr = STM\(rf_trans_real - r_arr_Mars)';
    
    r1_dot_new_trans = r1_dot - delta_r_dot_corr';
    r1_dot = r1_dot_new_trans;

    num_iter = num_iter + 1;
    
    end
    
    % Check convergence
    if num_iter>num_max_it
        warning('num MaxIter issue : The shooting method did not converge')
        r1_dot = [NaN NaN NaN];
        rf_dot_trans_real = [NaN NaN NaN];
        return
    end
end

[~, rf_dot_trans_real, ~, ~, ~] = FGKepler_dt(mu, r_dep_Earth, r1_dot, delta_t_target);

end