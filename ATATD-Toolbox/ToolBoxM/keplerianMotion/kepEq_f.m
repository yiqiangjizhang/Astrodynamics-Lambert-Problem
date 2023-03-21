function [f,errorflag] = kepEq_f(t,a,e,mu,f0,t0,option,imax,tol)

% kepEq_f.m - Keplerian Equation. It finds the true anomaly corresponding
%   to a certain time.
%
% PROTOTYPE:
%   [f,errorflag] = kepEq_f(t, a, e, mu, f0, t0, option, imax, tol)
%
% DESCRIPTION:
%   It finds the true anomaly corresponding to a time t. Multiple
%   revolutions are considered. Elliptic orbits (e < 1).
%
% INPUT:
%	t[1]            Time when the true anomaly is required [T].
%   a[1]            Semi-major axis [L].
%   e[1]            Eccentricity.
%   mu[1]           Planetary constant (mu = mass * G) [L^3/T^2].
%   f0[1]           Fixed true anomaly [rad].
%   t0[1]           Time corresponding to f0 [T]. Optional: default t0 = 0.
%   option[1]       Integer to select the printing options:
%                   1: a message is printed on the screen if the Newton
%                      loop does not converge within the required tolerance;
%                   0: no messages printed on the screen.
%                   > Optional: Default value = 0.
% 	imax[1]         Maximum number of iterations.
%                   > Optional: Default value = 5.
%  	tol[1]          Convergence tolerance for the Newton loop on
%                   abs(E-E0) [rad].
%                   > Optional: Default value = 1e-15.
%
% OUTPUT:
% 	f[1]            True anomaly at time t [rad] in [-Inf, Inf]. Note:
%                   Multiple revolutions are considered.
%	errorflag[1]    Integer identifying the error:
%                     0: No error;
%                     1: the Newton loop did not converge within the
%                        required tolerance;
%                     2: non elliptic orbit.
%
% CALLED FUNCTIONS:
%   (none)
%
% FUTURE DEVELOPMENT:
%   Extension to parabolic and hyperbolic orbits.
%
% AUTHOR:
%   Camilla Colombo, 28/11/2007, MATLAB, kepEq_f.m
%
% PREVIOUS VERSION:
%   Camilla Colombo, 25/11/2007, MATLAB, kep_eq.m
%       - multiple revolutions considered in the calculation of f.
%   Camilla Colombo, 28/11/2007, MATLAB, kep_eq1.m
%       - Header and function name in accordance with guidlines.
%
% CHANGELOG:
%   02/12/2008, Matteo Ceriotti: Modified Newton loop (in the case of a
%       limit cycle).
%   02/12/2008, REVISION: Matteo Ceriotti
%   30/12/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------

if nargin < 9
    tol = 1e-15;
    if nargin < 8
        imax = 5;
        if nargin < 7
            option = 0;
            if nargin < 6
                t0 = 0;
            end
        end
    end
end

% Eccentricity check
if e >= 1
    f = [];
    errorflag = 2;
    if option == 1
        fprintf('The orbit is not elliptic\n');
    end
    return
end

n  = sqrt(mu/a^3);  % [rad/T]

% Calculation of M0
% -pi < E0 < pi
E0 = acos((e+cos(f0))/(1+e*cos(f0)));
if sin(f0) < 0
    E0 = -E0;
end
M0 = E0-e*sin(E0);
% Test: if the keplerian parameters are from an asteroid, it is possible to
% check M0 by calculating M0_test
% [kep0_test,m0_test,M0_test] = NeoEphemeris(t0/86400,nNEO);

M = n*(t-t0)+M0;
nrev = fix((n*(t-t0))/(2*pi));

% Newton loop from Battin p. 217
E = M + e*sin(M)/(1-sin(M+e)+sin(M));   % Initial condition from Battin
                                        % p. 194, Eq. (5.4)
Mkold = NaN;
Mk = NaN;

niter = 0;
Eold = E + 2*tol;                       % To force at least one iteration
while abs(Eold-E) > tol && niter < imax
    niter = niter + 1;
    Eold = E;
    Mkoldold = Mkold;
    Mkold = Mk;
    Mk = E-e*sin(E);
    if Mk == Mkoldold % Newton loop is in a limit cycle
        % fprintf('Limit cycle. Re-setting Mk...\n');
        Mk = (Mk + Mkold)/2; % Set Mk in the middle of the limit cycle
    end
    E = E + (M-Mk)/(1-e*cos(E));
end

if abs(Eold-E) > tol
    errorflag = 1;
    if option == 1
        fprintf('The Newton loop did not converge within the required tolerance\n');
    end
else
    errorflag = 0;
end

f = 2*atan(sqrt((1+e)/(1-e)) * tan(E/2)); % [rad] [-pi pi]

if t >= t0 % Forward in time
    % f must follow f0
    if f < f0
        f = f + 2*pi;
    end
else % Backward in time
    % f must preceed f0
    if f > f0
        f = f - 2*pi;
    end
end

% Multiple revolutions
f = f + nrev*2*pi;

return