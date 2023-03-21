function out = rowNorm(in)

% rowNorm.m - Takes the 2-norm of each row.
%
% PROTOTYPE:
%   [out] = rowNorm(in)
%
% DESCRIPTION:
%   Takes the 2-norm of each row of the matrix in.
%
% INPUT:
%   in          Input matrix.
%
% OUTPUT:
%   out         2-norm, row by row.
%   
% CALLED FUNCTIONS:
%   (none)
%
% AUTHOR:
%   Matteo Ceriotti, 18/02/2010, MATLAB, rowNorm.m
%
% CHANGELOG:
%   18/02/2010, REVISION: Joan-Pau Sanchez
%
% -------------------------------------------------------------------------

out = sqrt(sum(abs(in).^2, 2));
