function hobj = plotOrbit(varargin)

% plotOrbit.m - Plot of the orbit.
%
% PROTOTYPE:
%   hobj = plotOrbit(kep5,mu)
%   hobj = plotOrbit(...,'PropertyName',PropertyValue,...)
%   hobj = plotOrbit(axes_handle,...)
%
% DESCRIPTION:
%   Plots a keplerian orbit without integrating the trajectory, but using
%   the Keplerian elements [a,e,i,Om,om].
%
% INPUT:
%	kep5[5]     Vector of Keplerian parameters of the orbit [a,e,i,Om,om]
%               a in [L], angles in [rad].
%   mu          Planetary gravity constant (mu = mass * G) [L^3/T^2].
%
% OUTPUT:
%   hobj        Handle to the plot object.
%
% CALLED FUNCTIONS:
%   kep2car
%
% AUTHOR:
%   Camilla Colombo, 17/02/2006, MATLAB, plotOrbit.m
%
% PREVIOUS VERSION:
%   Camilla Colombo, 17/02/2006, MATLAB, plotorbit.m
%       - Header and function name in accordance with guidelines.
%   
% CHANGELOG:
%   12/02/2007, Camilla Colombo: Changed kep2cart prototype
%   20/11/2007, REVISION: Matteo Ceriotti
%   20/11/2007, Matteo Ceriotti: - Added handle to plot object as output
%                                - Added initialisation of cart matrix for
%                                  speed
%   16/10/2009, Joan Pau Sanchez: Added the figure/axis handle.
%   02/11/2009, Joan Pau Sanchez: Added object property handling.
%   08/01/2010, Camilla Colombo: Header and function name in accordance
%       with guidelines.
%   05/10/2010, Camilla Colombo: Modified such that can work with color
%       string as single input.
%
%--------------------------------------------------------------------------

[inputs, properties]=parseparams(varargin);
% 1. Preraring Figure Object
if length(inputs)<3
    HAXIS = gca;
     kep5 = inputs{1};
       mu = inputs{2};
elseif ishandle(inputs{1})==0
        msg = ['The figure handle is not valid'];
        eid = sprintf('TOOLBOX:%s:propertyError', mfilename);
        error(eid,'%s',msg)
else
    try
        HAXIS=gca(inputs{1});
    catch
        HAXIS=inputs{1};  
    end
    hold on
    kep5 = inputs{2};
   	  mu = inputs{3};
end
hold on
%--------------------------------------------------------------------------
kep5 = kep5(1:5);

rdeg = pi/180;
cart = zeros(360,6);
for i = 1:360
    f = i*rdeg;
    cart(i,:) = kep2car([kep5 f],mu);
end

hobj = plot3(HAXIS,cart(:,1),cart(:,2),cart(:,3));   

xlabel('x');
ylabel('y');
zlabel('z');
hfig = gcf; 

if iscell(properties)
    if length(properties) == 1
        set(hobj,'Color',properties{1}) 
    else
        for i=1:2:length(properties)
            set(hobj,properties{i},properties{i+1}) 
        end
    end
end

return