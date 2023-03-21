function dy = pcr3bp_3D(t, x, mu,EphAn)

%        xdot = pcr3bp(t,x) ;
%
% with the LARGER MASS, m1 to the left of the origin at (-mu,0)
% and m2, or the planet (ie. Earth), is at (1 - mu, 0)
%
%                L4
% -L3----m1--+-----L1--m2--L2-
%                L5
%
% Shane Ross (revised 2.19.04)


mu1 = 1-mu; % mass of larger  primary (nearest origin on left)
mu2 =   mu; % mass of smaller primary (furthest from origin on right)

r3= ((x(1)+mu2)^2 + x(2)^2 + x(3)^2)^1.5;     % r: distance to m1, LARGER MASS
R3= ((x(1)-mu1)^2 + x(2)^2 + x(3)^2)^1.5;     % R: distance to m2, smaller mass

Ux = - x(1) + mu1*(x(1)+mu2)/r3 + mu2*(x(1)-mu1)/R3 ;
Uy = - x(2) + mu1* x(2)     /r3 + mu2* x(2)     /R3 ;
Uz = mu1* x(3)/r3 + mu2* x(3)/R3;

xdot = zeros(4,1);
xdot(1) = x(4);
xdot(2) = x(5);
xdot(3) = x(6);

xdot(4) = 2*x(5) - Ux ;
xdot(5) =-2*x(4) - Uy ;
xdot(6) = -Uz;

dy = xdot;

end
