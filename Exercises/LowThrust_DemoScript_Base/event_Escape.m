function [value,isterminal,direction]=event_Escape(t,x)
%--------------------------------------------------------------------------
muEarth = getAstroConstants('Earth','mu');
%--------------------------------------------------------------------------
% Current velocity
v=sqrt(x(4)^2+x(5)^2+x(6)^2);
% Escape velocity at altitude
vesc=sqrt(2*muEarth/sqrt(x(1)^2+x(2)^2+x(3)^2));

value = v-vesc; % when s/c reaches scape velocity value will be zero
isterminal = 1;
direction = 0;
end

                
