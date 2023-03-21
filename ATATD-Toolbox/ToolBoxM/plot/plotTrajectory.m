function [yend,y,t] = plotTrajectory(r0,v0,tof,mu,s)

% plotTrajectory.m - Plot of keplerian trajectory.
%
% PROTOTYPE:
%   [yend,y,t] = plotTrajectory(r0, v0, tof, mu[, s)
%
% DESCRIPTION:
%   Plots a keplerian trajectory of a body around a cental mass, given the
%   initial position and velocity in cartesian coordinates, and the final
%   time.
%
% Units of measure consistent each other.
%
% INPUT:
%   r0[1,3]     Cartesian initial position [L].
%   v0[1,3]     Cartesian initial velocity [L/T].
%   tof[1]      Time of flight [T].
%   mu[1]       Gravity constant of the central body [L^3/(M*T^2)].
%   s           Character string made from one element from the following
%               column to define the color of the line plotted
%                   b     blue
%                   g     green
%                   r     red
%                   c     cyan
%                   m     magenta
%                   y     yellow
%                   k     black
%               Optional. Default: s = 'b'
%
% OUTPUT:
%   yend[1,6]   Final state of the trajectory (position, velocity) [L, L/T].
%   y[n,6]     	State vector of the trajectory (position, velocity) [L, L/T].
%   t[n]        Time vector of the trajectory [T].
%
% CALLED FUNCTIONS:
%   dyn_2BP
%
% AUTHOR:
%   Matteo Ceriotti, 06/02/2007
%
% PREVIOUS VERSION:
%   Matteo Ceriotti, 06/02/2007, MATLAB, plot_trajectory.m
%       - Header and function name in accordance with guidlines.
%
% CHANGELOG:
%   06/02/2007, REVISION: Camilla Colombo
%   24/11/2007, Camilla Colombo: s added as an input, yend added as output
%   31/01/2009, Camilla Colombo: Substituted two_body_dynamics with d2b_eq.
%   05/10/2010, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------

if nargin < 5
    s = 'b';
end

% Initialises figure
line(0,0,0,'Marker','*','Color','k');

% Initialise integrator used to plot trajectories
ode45options = odeset('abstol',1e-11,'reltol',1e-9);
[t,y] = ode45(@dyn_2BP,[0, tof],[r0,v0],ode45options,mu);
yend = y(end,:);

line(y(:,1),y(:,2),y(:,3),'Color',s);
line(y(1,1),y(1,2),y(1,3),'Marker','x','Color',s) % Marker at departure point
line(yend(1),yend(2),yend(3),'Marker','o','Color',s) % Marker at arrival point


return
