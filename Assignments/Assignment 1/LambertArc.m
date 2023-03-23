function [r1_dot, rf_dot_trans_real] = LambertArc(mu, r_dep_Earth, r_arr_Mars, t_m, delta_t_target)

    % LambertArc.m - Lambert Arc Method for the Continuation Method of the Lambert Problem
    %
    % PROTOTYPE:
    %	[r1_dot, rf_dot_trans_real] = LambertArc(mu, r_dep_Earth, r_arr_Mars, t_m, delta_t_target)
    %
    % DESCRIPTION:
    %   This function is a Lambert Arc Method algorithm that uses the Continuation Method to solve the Lambert Problem.
    %   The Lambert Arc Method is a method to solve the Lambert Problem.
    %   The Continuation Method is a method to solve the Lambert Problem with a minimum energy transfer.
    %
    % INPUT:
    %   mu[1]               Gravitational parameter of the Sun [km3/s2]
    %   r_dep_Earth[3]      Position vector of the spacecraft at departure from the Earth [km]
    %   r_arr_Mars[3]       Position vector of the spacecraft at arrival to Mars [km]
    %   t_m[1]              Transfer method (+1, -1) for short and long path respectively
    %   delta_t_target[1]   Time of flight of the target [s]
    %
    % OUTPUT:
    %   r1_dot[3]            Velocity vector of the spacecraft at departure from the Earth [km/s]
    %   rf_dot_trans_real[3] Velocity vector of the spacecraft at arrival to Mars [km/s]
    %
    % CALLED FUNCTIONS:
    %   MinETransfer, FGKepler_dt, STM_Lambert
    %
    % AUTHOR:
    %   Yi Qiang Ji, 23/03/2023, MATLAB, LambertArc.m
    %
    % -------------------------------------------------------------------------

    % Minimum Energy Transfer of the spacecraft
    [~, ~, delta_t_min, r1_dot, ~] = MinETransfer(mu, r_dep_Earth, r_arr_Mars, t_m);

    % Evaluate tolerance
    tol = 1e-3; % [km]

    % Number of iterations
    num_it = 10; % Number of iterations until reach the desired slution

    % For loop of the Continuation Method
    for it = 1:num_it

        % Lambda parameter
        Lambda = it / num_it;
        delta_t = Lambda * delta_t_target + (1 - Lambda) * delta_t_min;

        % Number of maximum iterations
        num_max_it = 25;
        num_iter = 0; % Initialize number of iterations
        error = 1e19; % Initialize error

        % Loop
        while (error >= tol) && (num_iter < num_max_it)

            % Semi-major axis
            a = mu / (2 * mu / norm(r_dep_Earth) - norm(r1_dot) ^ 2);

            % Check sma (semi-major axis is negative)
            if a < 0
                %         warning('a seed error: The shooting method did not converge');
                r1_dot = [NaN NaN NaN];
                rf_dot_trans_real = [NaN NaN NaN];
                return
            end

            % Compute the real final position and F & G constants
            [rf_trans_real, ~, delta_E, F, G] = FGKepler_dt(mu, r_dep_Earth, r1_dot, delta_t);

            % Compute error
            error = norm(rf_trans_real - r_arr_Mars);

            % Compute the State Transition Matrix
            [STM] = STM_Lambert(mu, rf_trans_real, delta_E, F, G, r_dep_Earth, r1_dot, delta_t);
            delta_r_dot_corr = STM \ (rf_trans_real - r_arr_Mars)';

            % Apply correction
            r1_dot_new_trans = r1_dot - delta_r_dot_corr';
            r1_dot = r1_dot_new_trans;

            % Update number of iterations
            num_iter = num_iter + 1;

        end

        % Check convergence
        if num_iter > num_max_it
            %         warning('num MaxIter issue : The shooting method did not converge')
            r1_dot = [NaN NaN NaN];
            rf_dot_trans_real = [NaN NaN NaN];
            return
        end

    end

    % Compute the real final velocity
    [~, rf_dot_trans_real, ~, ~, ~] = FGKepler_dt(mu, r_dep_Earth, r1_dot, delta_t_target);

end
