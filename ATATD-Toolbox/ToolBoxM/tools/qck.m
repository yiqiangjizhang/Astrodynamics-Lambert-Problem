function angle = qck(angle)

% qck.m - Angle reduction within [0 2*pi] rad
%
% PROTOTYPE:
%   angle = qck(angle)
%
% DESCRIPTION:
%   Takes any angle and reduces it, if necessary, so that it lies in the
%   range from 0 to 2*pi radians.
%
% INPUT:
%   angle[1]    Ange to be reduced [rad].
%
% OUTPUT:
%   angle[1]    Angle reduced, if necessary, to the range from 0 to 2*pi
%               rad [rad].
%   
% CALLED FUNCTIONS:
%   (none)
%
% AUTHOR:
%   W.T. Fowler, July 1978, MATLAB, qck.m
%
% CHANGELOG:
%   20/08/1990, REVISION: Darrel Monroe
%   01/10/2007, REVISION: Matteo Ceriotti
%   30/12/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------

twopi = 2*pi;
 
diff = twopi * (fix(angle/twopi) + min([0,sign(angle)]));

angle = angle -diff;

return
