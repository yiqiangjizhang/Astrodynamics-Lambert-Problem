%% PlotEarth_Equatorial_frame.m
RE=getAstroConstants('Earth','radius');
% 2. Plot Axes
scale = 1;
figure1 = figure('Color',[1 1 1]);
% Create axes
axes1 = axes('Parent',figure1,'ZColor',[1 1 1],'YColor',[1 1 1],...
    'XColor',[1 1 1]);
hold(axes1,'all');
%--------------------------------------------------------------------------
[hfig] = drawPlanet('Earth',[0 0 0],figure1,scale);
% [hfig] = drawPlanet('Moon',[xEarth 0 0],figure1,scale*4);
%--------------------------------------------------------------------------
axis equal
%--------------------------------------------------------------------------
% Plot Axis
view(38,16)
%--------------------------------------------------------------------------
% X Axis
CMPos=[0 0 0];
p1=[5*RE 0 0];
%    p2=[xEarth-0.3 0 0];
line([ CMPos(1) p1(1)],[ CMPos(2) p1(2)],[ CMPos(3) p1(3)],'Color',[0 0 0])
%--------------------------------------------------------------------------
% Arrow head
% scale factor
H=gca;
xr=get(H,'xlim');
yr=get(H,'ylim');
zr=get(H,'zlim');
s=norm([xr(2)-xr(1),yr(2)-yr(1),zr(2)-zr(1)]);

r=s/300; % radius
h=6*r; % height
n=20; m=20; % grid spacing
[x,y,z]=cylinder(linspace(0,r,n),m);

H=surf(x+p1(1),y+p1(2),z*h+p1(3));
set(H,'FaceColor','k','EdgeColor','k')
P=p1-CMPos; % position vector
U=P/norm(P); % unit vector
D=cross([0,0,1],-U); % rotation axis
if norm(D)<eps, D=[1,0,0]; end
A=acos(dot([0,0,1],-U)); % rotation angle
rotate(H,D,A*180/pi,p1)
%--------------------------------------------------------------------------
% Y axis
CMPos=[0 0 0];
py=[0 5*RE 0];
%    p2=[xEarth-0.3 0 0];
line([ CMPos(1) py(1)],[ CMPos(2) py(2)],[ CMPos(3) py(3)],'Color',[0 0 0])
% Arrow head Y
% scale factor
H=gca;
xr=get(H,'xlim');
yr=get(H,'ylim');
zr=get(H,'zlim');
s=norm([xr(2)-xr(1),yr(2)-yr(1),zr(2)-zr(1)]);

r=s/300; % radius
h=6*r; % height
n=20; m=20; % grid spacing
[x,y,z]=cylinder(linspace(0,r,n),m);

H=surf(x+py(1),y+py(2),z*h+py(3));
set(H,'FaceColor','k','EdgeColor','k')
P=py-CMPos; % position vector
U=P/norm(P); % unit vector
D=cross([0,0,1],-U); % rotation axis
if norm(D)<eps, D=[1,0,0]; end
A=acos(dot([0,0,1],-U)); % rotation angle
rotate(H,D,A*180/pi,py)
%--------------------------------------------------------------------------
% Z axis
CMPos=[0 0 0];
pz=[0 0  5*RE];
%    p2=[xEarth-0.3 0 0];
line([ CMPos(1) pz(1)],[ CMPos(2) pz(2)],[ CMPos(3) pz(3)],'Color',[0 0 0])
% Arrow head Y
% scale factor
H=gca;
xr=get(H,'xlim');
yr=get(H,'ylim');
zr=get(H,'zlim');
s=norm([xr(2)-xr(1),yr(2)-yr(1),zr(2)-zr(1)]);

r=s/300; % radius
h=6*r; % height
n=20; m=20; % grid spacing
[x,y,z]=cylinder(linspace(0,r,n),m);

H=surf(x+pz(1),y+pz(2),z*h+pz(3));
set(H,'FaceColor','k','EdgeColor','k')
P=pz-CMPos; % position vector
U=P/norm(P); % unit vector
D=cross([0,0,1],-U); % rotation axis
if norm(D)<eps, D=[1,0,0]; end
A=acos(dot([0,0,1],-U)); % rotation angle
rotate(H,D,A*180/pi,pz)
%--------------------------------------------------------------------------