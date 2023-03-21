function [rp,rho_bpl] = rho2rp(rho_car,vA,vB,muA)

% rho2rp.m - Conversion from relative position to pericentre of the
%   hyperbola
%
% PROTOTYPE:
%   rp = rho2rp(rho_car,vA,vB,mu)
%
% DESCRIPTION:
%	Computes the future/past pericentre distance rp of body B from body A
%   given the relative distance rho_car (in Cartesian), the absolute
%   velocity vector of body A, vA (in Cartesian), and the absolute velocity
%   vector of body B, vB (in Cartesian), and the gravity constant of body A.
%   Note: When using this function for asteroid deflection, we usually use
%       the nominal unperturbed velocity of body B but for a more precise
%       result, the velocity on the deflected orbit of B should be used
%       instead.
%
% INPUT:
%   rho_car     Relative position between A and B (Cartesian).
%  	vA          Absolute velocity vector of A (Cartesian).
%   vB          Absolute velocity vector of B (Cartesian).
%   muA         Gravity constant of body A.
%
%  OUTPUT:
% 	rp          Pericentre distance on the hyperbola centered on body A.
%   rho_bpl     Relative position between A and B expressed in the b-plane
%               reference system. The b-parameter vector on the b-plane is:
%               b = [rho_bpl(1),rho_bpl(3)]
%
% CALLED FUNCTIONS:
%   car2bpl
%
% EXAMPLE: Body A Earth, body B asteroid. rho_car can be the relative
%   position of the deflected asteroid with respect to its nominal position
%   on the undeflected orbit plus the MOID distance: dr+Dr.
%
% AUTHOR:
%   Massimiliano Vasile, 20/11/2009
%
% CHANGELOG:
%   08/11/2009, REVISION, Camilla Colombo & Joan Pau Sanchez.
%   18/01/2010, Camilla Colombo: The relative position between A and B has
%       to be inputed in Cartesian reference frame.
%   03/06/2010, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%   18/11/2010, Camilla Colombo: Added rho_bpl as output.
%
% -------------------------------------------------------------------------


U_car = vB-vA;	% Relative velocity of body B with respect to body A.

rho_bpl = car2bpl(rho_car,U_car,vA); % Cartesian to b-plane reference frame.

rp = sqrt(muA^2/norm(U_car)^4+norm(rho_bpl(1:2:3))^2)-muA/norm(U_car)^2;

% rho_bpl(1:2:3) takes the component 1 and 3 of rho_bpl:
% b = [rho_bpl(1),rho_bpl(3)]

return
