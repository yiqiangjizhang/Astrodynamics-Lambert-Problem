function x_bpl = car2bpl(x_car,U_car,vp_car)

% car2bpl.m - Vector reference frame transformation.
%	Cartesian to b-plane reference frame.
%
% PROTOTYPE:
%	x_bpl = car2bpl(x_car, U_car, vp_car)
%
% DESCRIPTION:
%   Transformation from Cartesian to b-plane reference frame.
%   Cartesian (car) reference frame: {x,y,z}
%       inertial reference frame.
%   b-plane reference frame: {xi,eta,zeta}
%   The b-plane is the plane perpendicular to the incoming relative
%   velocity of the small body at the planet arrival, containing the planet.
%       eta-axis: direction of the incoming relative velocity of the small
%                 body on its arrival.
%       zeta-axis: in the b-plane, direction opposite to the projection of
%                  the heliocentric velocity of the planet on the b-plane.
%       xi-axis: in the b-plane, completes the reference frame.
%
% INPUT:
%	x_car[3]    Vector to be transformed, expressed in {x,y,z}.
%  	U_car[3]    Velocity [L/T] of the small body relative to the planet,
%               expressed in {x,y,z}.
% 	vp_car[3]   Orbital velocity of the planet [L/T], expressed in {x,y,z}.
%
% OUTPUT:
% 	x_bpl[3,1]	Vector transformed in {xi,eta,zeta}.
%
% CALLED FUNCTIONS:
%   crossFast
%
% AUTHOR:
%   Camilla Colombo, 04/05/2007, MATLAB, car2bpl.m
%
% PREVIOUS VERSION:
%   Camilla Colombo, 04/05/2007, MATLAB, car_bplT.m
%       - Header and function name in accordance with guidlines.
%   
% CHANGELOG:
%   21/05/2007, REVISION: Matteo Ceriotti
%   11/02/2008, Matteo Ceriotti: Help improved.
%   30/09/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%   29/03/2010, Camilla Colombo: crossFast used.
%
% -------------------------------------------------------------------------

x_car  = x_car(:);
U_car  = U_car(:);
vp_car = vp_car(:);

nn = U_car/norm(U_car);
ee = crossFast(vp_car,nn)/norm(crossFast(vp_car,nn));
cc = crossFast(ee,nn);

x_bpl = [ee nn cc]'*x_car;

return