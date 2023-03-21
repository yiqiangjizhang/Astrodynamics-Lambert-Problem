function semip = semiperiodHohm(mu,a)

% semiperiodHohm.m - Semi period of an Hohmann transfer.
%
% PROTOTYPE:
%   semip = semiperiodHohm(mu,a)
%
% DESCRIPTION:
%   Calculates the semi period of an Hohmann transfer given the semi-major
%   axis.
%
% INPUT:
%	mu      Planetary constant (mu = mass * G) [L^3/T^2].
%	a       Semi-major axis [L].
%
% OUTPUT:
%	semip	Semi-period of the orbit [T].
%
% CALLED FUNCTIONS:
%   (none)
%
% AUTHOR:
%   Camilla Colombo, 02/08/2006, MATLAB, semiperiodHohm.m
%
% PREVIOUS VERSION:
%   Camilla Colombo, 02/08/2006, MATLAB, semiperiod_hohm.m
%       - Header and function name in accordance with guidlines.
%   
% CHANGELOG:
%   14/02/2007, Camilla Colombo: name changed (THohm1.m)
%   04/11/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------


semip = pi*sqrt(a^3/mu);
return