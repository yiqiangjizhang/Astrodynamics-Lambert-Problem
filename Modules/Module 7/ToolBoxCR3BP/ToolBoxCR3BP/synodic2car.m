function [X_c, V_c] = synodic2car(x, v, time, mu, theta_0)

    %Converts Synodic RF to Cartesian Inertial
    %Primary Sun, Secondary Earth (ephSS_car)
    
    %addpath(genpath('C:\Users\Rita\Dropbox\PHD\Matlab\CranRepo'))
    
    if isempty(theta_0)
    	theta_0 = 0;
    end
    
    %Rotation
	theta=time - theta_0;
    
    A_t = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

    %Renormalize and change center
    v = (v + [0 mu 0]);                                                  
    x = (x + [mu 0 0]);  
    
    V_c = (A_t*[v(1)-x(2) v(2)+x(1) v(3)]')';
    X_c = (A_t*x')';
       
end