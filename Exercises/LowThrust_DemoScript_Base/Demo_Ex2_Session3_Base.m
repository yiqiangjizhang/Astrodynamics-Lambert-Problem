clear all; close all; clc; 
%% Demo_Ex2_Session3.m
% Propagate, using Matlab ODE solvers, the escape trajectory described in
% the previous exercise using a constant thrust in tangential direction.
% Initial Orbit: Circular 250 km altitude parking orbit. Apogee motor
% specs: Thrust 425 N,  Isp 320s. Spacecraft wet mass : 2000 kg.
%% Compute how long you need to thrust for to reach escape velocity?
% Note: use equations in slide 10.
% DT?



%% Low Thrust Propagation
% Use the above computed thrusting time DT, as propagation time for the ODE

% Spacecraft in circular orbit 250 km
% RE - radius of the Earth
% Vc - velocity circular orbit
% Msc- spacecraft mass

RE = 6378; % [km]
h = 250; % [km]
mu_E = 
Vc = sqrt(mu_E/(RE + h)); 
Msc = 2000; % [kg]

SV0=[RE+250 0 0 0 Vc 0 Msc];

% LOW THRUST Controls
        % Tangential Thrust Vector (Same direction as velocity)
        Tx=@(t,x) (Tmax*x(4)/sqrt(x(4)^2+x(5)^2+x(6)^2));
        Ty=@(t,x) (Tmax*x(5)/sqrt(x(4)^2+x(5)^2+x(6)^2));
        Tz=@(t,x) (Tmax*x(6)/sqrt(x(4)^2+x(5)^2+x(6)^2)); %T=Tmx*V/|V|

        %Equation of motion ¨r+mu*r/R^3 = T/m
        F=@(t,x)   [x(4); %dx/dt=Vx
                    x(5); %dy/dt=Vy
                    x(6); %dz/dt=Vz
                    -muEarth*x(1)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3)+Tx(t,x)/(x(7)*1000); %dVx/dt
                    -muEarth*x(2)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3)+Ty(t,x)/(x(7)*1000); %dVx/dt
                    -muEarth*x(3)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3)+Tz(t,x)/(x(7)*1000); %dVx/dt
                    -sqrt(Tx(t,x)^2 + Ty(t,x)^2 + Tz(t,x)^2)/ve]; %dm/dt=T/Isp;

        options=odeset('RelTol',1e-6,'AbsTol',1e-7,'Refine',50);
        
        %Propagation
        [T,X]=ode45(F,[0,DT],SV0,options);

%--------------------------------------------------------------------------
% Plot reference frame
run PlotEarth_Equatorial_frame.m
plot3(X(:,1),X(:,2),X(:,3))
view(0,90)

%% Did we reach escape velocity?
% A spacecraft will be in a escape trajectory if at any given altitude it
% has a velocity equivalent to sqrt(2*muEarth/r);
% Hence, let's check the final velocity of the above propagation.



%% Propagate the trajectory until reaching escape conditions
% Use MATLAB’s ODE event feature to stop the propagation when it reaches
% escape velocity. For information on how to do this go to:
% https://uk.mathworks.com/help/matlab/math/ode-event-location.html


% Spacecraft in circular orbit 250 km
% RE - radius of the Earth
% Vc - velocity circular orbit
% Msc- spacecraft mass
SV0=[RE+250 0 0 0 Vc 0 Msc];
tf=10*DT; % make sure you propagate for substantially longer than the 
          % previous propagation time

% LOW THRUST Controls
        % Tangential Thrust Vector (Same direction as velocity)
        Tx=@(t,x) (Tmax*x(4)/sqrt(x(4)^2+x(5)^2+x(6)^2));
        Ty=@(t,x) (Tmax*x(5)/sqrt(x(4)^2+x(5)^2+x(6)^2));
        Tz=@(t,x) (Tmax*x(6)/sqrt(x(4)^2+x(5)^2+x(6)^2)); %T=Tmx*V/|V|

        %Equation of motion ¨r+mu*r/R^3 = T/m
        F=@(t,x)   [x(4); %dx/dt=Vx
                    x(5); %dy/dt=Vy
                    x(6); %dz/dt=Vz
                    -muEarth*x(1)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3)+Tx(t,x)/(x(7)*1000); %dVx/dt
                    -muEarth*x(2)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3)+Ty(t,x)/(x(7)*1000); %dVx/dt
                    -muEarth*x(3)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3)+Tz(t,x)/(x(7)*1000); %dVx/dt
                    -sqrt(Tx(t,x)^2 + Ty(t,x)^2 + Tz(t,x)^2)/ve]; %dm/dt=T/Isp;
             
        
        options=odeset('RelTol',1e-6,'AbsTol',1e-7,'Refine',50,'Events',@(t,x)event_Escape(t,x));
        % event_Escape is defined as in slide 21
        
        %Propagation
        [T,X,TE,YE,IE]=ode45(F,[0,tf],SV0,options);

%--------------------------------------------------------------------------
% Plot reference frame
run PlotEarth_Equatorial_frame.m
plot3(X(:,1),X(:,2),X(:,3))
view(0,90)

%--------------------------------------------------------------------------
% Is the final state vector in a escape trajectory? 

%--------------------------------------------------------------------------
% What is the gravity loss that we will have if we intent to thrust until
% escape velocity from parking orbit with our engine?

% the total thrust time is TE

% Again, use the equations in slide 10.
