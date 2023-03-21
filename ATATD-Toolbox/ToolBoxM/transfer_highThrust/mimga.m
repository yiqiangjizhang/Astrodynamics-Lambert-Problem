function [dv,vin,error,xp2,vp2,output] = mimga(t0, T, s, multirev, v0, options ,sys)

% mimga.m - Multi gravity assists.
%
% PROTOTYPE:
%   [dv,vin,error,xp2,vp2,output] = mimga(t0, T, s, multirev, v0, options, sys)
%
% DESCRIPTION:
%   It computes the dv for a sequence of gravity assists and deep space
%   manoeuvres, in order to reach a given senquence of planets.
%   It uses the traditional patched-conic approximation, and the main body
%   during deep space flights is the Sun. It uses a multi revolution
%   lambert solver. Only flybys of planets (1 <= s <= 9) are allowed.
%
% INPUT:
%   t0[1]       Departure time [d, MJD2000].
%   T           Vector of parameters.
%               If options(1)=0:
%                   T = [alpha1, TOF1,
%                       gamma2, rp2, alpha2, TOF2,
%                       gamma3, rp3, alpha3, TOF3,
%                       ...];
%               If options(1)=1:
%                   T = [gamma1, rp1, alpha1, TOF1,
%                       gamma2, rp2, alpha2, TOF2,
%                       ...];
%               If options(2)=0:
%                   T = [...,
%                       gamman, rpn, alphan, TOFn];
%               If options(2)=1:
%                   T = [...,
%                       gamman, rpn, alphan, TOFn,
%                       gamma, rp];
%               Where:
%                   TOF: Time of flight [d];
%                   gamma: angle of hyperbola plane [rad];
%                   rp: radius of the pericentre of the hyperbola
%                       (adimensionalised with the radius of the planet).
%                   alpha: fraction of TOF before the deep space manoeuvre.
%   v0[1,3]     Initial velocity vector [km/s]. It can be relative to the
%               first planet or absolute, depending on options(3).
%   s           Planetary sequence. The first planet is the departure
%               planet, the last planet is the arrival planet. The flyby of
%               the starting and the ending planet is regulated by
%               options(1:2). All other planets will be flown by. Depending
%               of the options, s must contain at least one or two planets.
%   multirev    Matrix containing, in each column (relative to each Lambert
%               arc), the type of transfer (direct or retrograde), the
%               number of revolutions for the Lambert arc, the case of
%               transfer (low or high energy):
%                   multirev = [type1, type2,  ...;
%                               nrev1, nrev2, ...;
%                               case1, case2, ...];
%               Leave it empty [] for using 0 rev, low energy, direct
%               transfer for each arc. Or it's possible to fill only the
%               first 1 or 2 rows, to take 0 rev and low energy.
%   options     Vector of options. It is optional, and it is also possible
%               to specify only the first options on the left of the vector.
%               All the others are assumed 0.
%               options(1) determines how to start the journey.
%                   0:  The journey starts at s(1) with a deep space flight,
%                       with departure velocity v0. No swing-by is
%                       performed at body s(1).
%                   1:  The journey starts with a swing-by at planet s(1).
%                       v0 is the incoming velocity before the swing-by
%                       manoeuvre.
%               options(2) determines how to end the journey.
%                   0:  The journey ends at the encounter of the last planet,
%                       s(end) without doing the swing-by of it.
%                   1:  The function ends after a swing-by of the last
%                       planet, s(end).
%               options(3) is used to know whether v0 is absolute or
%                       relative to the first planet of the sequence.
%                   0:  v0 is the absolute initial velocity of the spacecraft,
%                       given in cartesian coordinates.
%                   1:  v0 input is relative to the first planet, in
%                       cartesian coordinates.
%                       This case can be used with v0 = [0 0 0], and in
%                       conjunction with alpha1=0, to create a trajectory
%                       in which there is no real DSM in the first leg. In
%                       fact, the DSM is moved to emulate launch. Note that
%                       in this case, the first element of the output dv 
%                       will still be considered as a DSM.
%                   2:  v0 input is relative to the first planet, and given
%                       as the modulus, the in-plane angle with respect to
%                       the planetary velocity, and the out-of-plane angle
%                       [rad].
%                   3:  As in case 2, but the 2 angles are given adimensionally
%                       in the interval [0, 1], such that the sphere is sampled
%                       uniformly. See Weisstein, Eric W. "Sphere Point
%                       Picking." From MathWorld--A Wolfram Web Resource.
%                       http://mathworld.wolfram.com/SpherePointPicking.html
%               options(4) is used to plot the trajectory or print data.
%                   0: No plot.
%                   1: The trajectory is plotted in a new figure.
%                   2: Plot and print trajectory data on the screen.
%               options(5) is used to ask the output structure containing data
%                   about the trajectory. Since generating this structure
%                   requires further calculations, it's reccomended not to ask
%                   if not needed. Anyway, the data are not computed if output
%                   is not in the output list.
%                   1: output structure generated (if in the output list).
%                   0: output structure not generated.
%   sys         Optional. A structure which defines the planetary system and
%               ephemerides to use:
%                   .mu: Gravitational constant of the main attrector;
%                   .mu_planets: Gravitational constants of the satellites;
%                   .R_planets: Radii of the satellites;
%                   .eph_handle: Handle to the ephemeris function for the
%                       bodies, with prototype: [x, v] = eph_handle(id, t).
%                   .eph_kep_handle: Handle to the ephemeris function for the
%                       keplerian parameters of the bodies, with prototype:
%                       kep = eph_kep_handle(id, t).
%               If not specified or empty, uses the Solar System
%               (astroConstants and EphSS_car).
%
% OUTPUT:
%   dv          Vector of DSM dv [km/s]. No launch and no brake dv is
%               computed.
%   vout        Absolute velocity vector at the end of the sequence [km/s].
%   error       0 if and only if everything ok. If error occurs, dv = 1e6.
%   xp2         Vector of the position of the last encountered planet [km].
%   vp2         Vector of the velocity of the last encountered planet [km/s].
%   output      Output structure with the trajectory data. Generated only if
%               options(5)=1.
%
% EXAMPLES:
%   Launch from Earth and deep space flight to Venus:
%       t0 = [3359.38387247647];
%       T = [0 170.49933001840];
%       s = [3,2];
%       [dv,vout,error,xp2,vp2] = mimga(t0, T,s,[],[0,0,0],[0,0,1]);
%   
%   Just a flyby of Venus:
%       t0 = [3359.38387247647+170.49933001840];
%       T = [3.26576512201 2.07228045370];
%       s = [2];
%       vin = [... ... ...];
%       [dv,vout,error,xp2,vp2] = mimga(t0, T,s,[],vin,[1,0,0]);
%
%   Launch from Earth, deep space flight and flyby of Venus:
%       t0 = [3359.38387247647];
%       T = [0 170.49933001840 3.26576512201 2.07228045370];
%       s = [3,2];
%       [dv,vout,error,xp2,vp2] = mimga(t0, T,s,[],[0,0,0],[0,1,1]);
%
%   Complete EVEEJ transfer, with plot:
%       t0 = [3359.38387247647];
%       T = [0, 170.49933001840,...
%            3.26576512201, 2.07228045370, 0.06197382751, 320.21104182897, ...
%            2.65131587336, 1.66668311162, 0.07728172747, 730.48301546338, ...
%            3.16200529828, 1.41291334360, 0.10740729409, 847.09714061671];
%       s = [3,2,3,3,5];
%       [dv,vout,error,xp2,vp2] = mimga(t0,T,s,[],[0,0,0],[0,0,1,1]);
%
% CALLED FUNCTIONS:
%   lambertMR, swingby, kepPro, astroConstants, dyn_2BP, car2kep.
%   If sys is specified, ephemerides functions pointed by sys.eph and
%   sys.eph_kep; otherwise, EphSS_car and astroConstants.
%
% ORIGINAL VERSION:
%   Massimiliano Vasile, 2002, MATLAB, mimga.m
% 
% AUTHOR:
% Matteo Ceriotti, 11-11-2006, MATLAB, mimga.m (Freely inspired by
%   Massimiliano Vasile's mimga).
%
% PREVIOUS VERSION:
%   Matteo Ceriotti, 11-11-2006, MATLAB, mimga3_3.m
%       - Header and function name in accordance with guidlines.
%
% CHANGELOG:
%   01-08-2006, Matteo Ceriotti: Added support for inital velocity given
%       with modulus and 2 angles.
%   09-11-2006, Matteo Ceriotti: Does not read the radius and the
%     	mu of the planet anymore, if no flyby will be performed.
% 	10-11-2006, Matteo Ceriotti: When plotting, chooses the color
%     	as a function of the semimajor axis of the planet
%   Version 3, Matteo Ceriotti, 01-01-2007: Multi rev. Lambert
% 	13-02-2007, Matteo Ceriotti: Uses keppro3, cart2kep.
% 	23-07-2007, Matteo Ceriotti: When error occurs, dv = 1e6.
%   26-07-2007, Matteo Ceriotti: Added verbose possibility in options(4).
%   30-07-2007, Matteo Ceriotti: Added the possibility of uniform sampling
%       of the sphere for v0.
%   Version 3.3, Matteo Ceriotti, 05-06-2008: Added possibility of other
%    	planetary systems. Changed planet orbit colours in plotting.
%      	Modified way of plotting planets' orbits: now they are not
%      	integrated anymore, but ephemerides are called at each instant
%   	of time.
% 	01-08-2008, Matteo Ceriotti: Added the output structure as an output.
% 	19-12-2008, Matteo Ceriotti: Changed two_body_dynamics to d2b_eq.
% 	23/02/2009, Matteo Ceriotti: Fixed empty swingby field in
%     	output structure when no swing-by at beginning is performed.
%      	Added rp in swingby in output structure. Added body and relative
%       velocity in output structure, leg(1).part(1) when alpha(1)==0.
%       Added dimensional rp and hp.
%   04/10/2010, Camilla Colombo: Header and function name in accordance
%       with guidlines (also changed name of called functions: keppro3
%       to kepPro, astro_constants to astroConstants, d2b_eq to dyn_2BP,
%       cart2kep to car2kep)
%
% -------------------------------------------------------------------------

if nargin < 7
    sys = [];
    if nargin < 6
        options=[];
    end
end

% Completes options vector
if length(options)<5
    options = [options,zeros(1,5-length(options))];
end

% Defines the default planetary system (Solar System)
if isempty(sys)
    sys.mu = astroConstants(4); % Planetary constant of the Sun
    sys.mu_planets = astroConstants(11:19);  % Default: Solar System
    sys.R_planets = astroConstants(21:29); % Default: Solar System
    sys.eph_handle = @EphSS_car; % Default: Solar System
    sys.eph_kep_handle = @EphSS_kep; % Default: Solar System
end

if options(4)==2 % Verbose
    fprintf('  *** mimga3_3 - Verbose mode ***\n');
    fprintf('  ***\nInput parameters\n');
    fprintf('t0 = %f\nT =\n[',t0);
    for i=1:length(T)
        fprintf('%f ',T(i));
    end
    fprintf(']\ns =\n[ ');
    for i=1:length(s)
        fprintf('%d ',s(i));
    end
    fprintf(']\nmultirev = [\n');
    for j=1:size(multirev,1)
        for i=1:size(multirev,2)
            fprintf('%d ',multirev(j,i));
        end
        fprintf(';\n');
    end
    fprintf('v0 =\n[');
    for i=1:length(v0)
        fprintf('%f ',v0(i));
    end
    fprintf(']\noptions =\n[');
    for i=1:length(options)
        fprintf('%d ',options(i));
    end
    fprintf(']\n');
end

if options(4) % Plot flag
    % Initialises figure
    hold on;
    line(0,0,0,'Marker','*');
    xlabel('x, km'); ylabel('y, km'); zlabel('z, km');
    % *** AU=astroConstants(2);
    % Initialise integrator used to plot trajectories
    ode45options = odeset('abstol',1e-9,'reltol',1e-7);
end
% Completes multirev matrix
nlambertarcs=length(s)-1;
if size(multirev,1)==0 % Not specified anything
    multirev(1,:) = zeros(1,nlambertarcs); % Assume direct (1) transfer for each arc
end
if size(multirev,1)==1 % Not specified the number of revolutions
    multirev(2,:) = zeros(1,nlambertarcs); % Assume 0 revolutions for each arc
end
if size(multirev,1)==2 % Not specified the case
    multirev(3,:) = zeros(1,nlambertarcs); % Assume low energy case (0) transfer for each arc
end

error=0;
% *** mu=astroConstants(4); % Planetary constant of the Sun

planet1=s(1);
time=t0;
Tindex=1;
sindex=2;
[xp1, vp1] = sys.eph_handle(planet1, time); % Ephemeris of planet1

if options(4) % Plot flag
    line(xp1(1),xp1(2),xp1(3),'Marker','*') % Marker at departure point
end

if options(5) && nargout > 4 % Initialises output
    output.sys = sys;
    output.s = s;
    if options(1) == 1 % Initial swingby
        output.composition = 's';
    else
        output.composition = [];
    end
    for i = 1:length(s)-2
        output.composition = [output.composition, 'ds'];
    end
    if options(2) == 1 % Final swingby
        output.composition = [output.composition, 'ds'];
    else
        if length(s) > 1
            output.composition = [output.composition, 'd'];
        end
    end
else
    output = [];
end

dv = [];
switch options(3)
    case 0, % v0 is absolute, so nothing to do.

    case 1, % v0 is relative to the first planet, and cartesian
        v0=v0+vp1; % Conversion to absolute initial velocity
    case {2, 3}, % v0 is relative to the first planet, and given as modulus and 2 angles: the in-plane and the out-of plane.
        if options(3)==3 % A pre-transformation is needed, as the angles are given in adimensional components to uniformly sample the sphere
            % Uniform sampling of the sphere
            v0(2)=v0(2)*2*pi-pi/2;
            v0(3)=acos(2*v0(3)-1)-pi/2;
        end
        mod_v=v0(1); % Modulus of initial velocity relative to the departure planet
        psi=v0(2)-pi/2; % In-plane rotation
        theta=v0(3); % Out-of-plane rotation
        % Transposed Euler angles matrix in the special case phi=0.
        R_T=[cos(psi), -cos(theta)*sin(psi),  sin(theta)*sin(psi);...
             sin(psi),  cos(theta)*cos(psi), -sin(theta)*cos(psi);...
                0    ,       sin(theta)    ,       cos(theta)    ];
        v_rel_tnh=R_T*[0;mod_v;0]; % Relative velocity in the tangent-normal-binormal (h) reference frame
        v0=tnh2car(v_rel_tnh,[xp1,vp1])';
        v0=v0+vp1; % Conversion to absolute initial velocity
    otherwise
        disp('Error: bad options(3)!');
        return
end

vp2=vp1; % Just as initialisation. If there is only one planet, this is useful for the calculation of the braking dv.
xp2=xp1; % Just as initialisation. If there is only one planet, this is useful for the output ephemeris.

switch options(1)
    case 0,
        % Same as (2), but the v0 velocity vector is after the flyby.
        % Basically, the flyby has been already done, and so starts the deep
        % space flight.
        vout=v0; % v0 is the velocity at the end of the flyby
    case 1,
        % Using the velocity vector v0 at the planet encounter, starts a swing
        % by of the first planet in the list, then continues with deep space
        % manoeuvre, swingby on the second planet, and so on (till the last
        % planet).
        vin=v0; % v0 is the velocity before the first flyby
    otherwise
        disp('Error: bad options(1)!');
        return
end

% Remeber to initialise Tindex, sindex, time, vin, planet1, xp1, vp1, nflybys, dv before the for

% Flybys and deep space manoeuvres
while Tindex<=length(T)
    % Here planet1 is the flyby planet, and planet2 is the planet towards
    % which the spacecraft is pointing (if available).
    
    % Flyby of planet1
    if sindex==2 && options(1)==0 % sindex==2 means that we are at the first planet
        % Skip the flyby of the first planet.
    else
        % Flyby of planet1
        % Reads mu and radius of the planet to be flown by.
        mup = sys.mu_planets(planet1); % Gravity constant of the planet1
        rplanet = sys.R_planets(planet1); % Radius constant of the planet1

        % Reads flyby variables from vector T and s
        rp_nondim = T(Tindex+1); % Pericentre of the hyperbola (non-dimensional)
        rp=T(Tindex+1)*rplanet; % Pericentre of the hyperbola (dimensional)
        psi=T(Tindex+0); % Hyperbola plane angle
        Tindex=Tindex+2;
        
        % Swing-by
        vinf1=vin-vp1; % Velocity relative to the planet
        n=cross(vinf1,vp1); % Reference vector
        n=n/norm(n);
        vinf2=swingby(vinf1,rp,mup,psi,n); vinf2 = vinf2(:)'; % Swingby
        vout=vinf2+vp1; % Absolute velocity after the swingby
        if options(4)==2 % Verbose
            % Computes the deflection angle for verbosing
            theta_inf_verbose = acos(-mup/rp/(norm(vinf1)^2+mup/rp)); % See Kaplan "Modern spacecraft dynamics and control", pag. 93
            delta_verbose = 2*theta_inf_verbose - pi; % Deflection angle
            temp = car2rth(vinf1,[xp1,vp1]);
            fprintf('  ***\nSwingby of %d - time = %f d,MJD2K\n', planet1, time);
            fprintf('Relative incoming velocity in r-th-h (also called R-S-T) [km/s]:\n [%f %f %f]\n',temp(1),temp(2),temp(3));
            fprintf('Relative incoming velocity in x-y-z [km/s]:\n [%f %f %f]\n', vinf1(1),vinf1(2), vinf1(3));
            fprintf('|vinf| = %f km/s; Deflection angle (delta) = %f deg\n', norm(vinf1), delta_verbose*180/pi);
            fprintf('gamma = %f rad; rp = %f km (hp = %f km)\n',psi,rp,rp-rplanet);
            temp = car2rth(vinf2,[xp1,vp1]);
            fprintf('Relative outgoing velocity in r-th-h (also called R-S-T) [km/s]:\n [%f %f %f]\n',temp(1),temp(2),temp(3));
            fprintf('Relative outgoing velocity in x-y-z [km/s]:\n [%f %f %f]\n',vinf2(1),vinf2(2),vinf2(3));
        end
        if options(5) && nargout > 4 % Output requested
            if options(1) == 1 % Initial swing-by
                indexswingby = sindex - 1;
            else % No initial swing-by
                indexswingby = sindex - 2;
            end
            
            output.swingby(indexswingby).incoming.v_rel_rth = car_rthT(vinf2,[xp1,vp1]);
            output.swingby(indexswingby).incoming.v_rel_xyz = vinf1;
            output.swingby(indexswingby).incoming.kep = car2kep([xp1, vin], sys.mu);
            output.swingby(indexswingby).incoming.v = vin;
            output.swingby(indexswingby).outgoing.v_rel_rth = car_rthT(vinf2,[xp1,vp1]);
            output.swingby(indexswingby).outgoing.v_rel_xyz = vinf2;
            output.swingby(indexswingby).outgoing.kep = car2kep([xp1, vout], sys.mu);
            output.swingby(indexswingby).outgoing.v = vout;
            output.swingby(indexswingby).vinf = norm(vinf1);
            output.swingby(indexswingby).t = time;
            output.swingby(indexswingby).body = planet1;
            output.swingby(indexswingby).x = xp1; % Position of the planet (and the spacecraft in the heliocentric frame)
            output.swingby(indexswingby).vp = vp1; % Heliocentric velocity of the planet
            output.swingby(indexswingby).rp_nondim = rp_nondim; % Non dimensional radius of pericentre
            output.swingby(indexswingby).rp = rp; % Dimensional radius of pericentre [km]
            output.swingby(indexswingby).hp = rp - rplanet; % Altitude of closest approach [km]
        end
    end
    
    % Deep space flight (from planet1 to planet2) and manoeuvre
    if sindex>length(s) % sindex>length(s) means that we are at the last planet
        % Do nothing: only the flyby of the last planet had to be done.
        vin=vout;
    else
        % Reads trajectory variables from vector T and s
        alpha=T(Tindex+0); % Time fraction before deep space manoeuvre
        tof=T(Tindex+1); % Time of flight to the next planet
        planet2=s(sindex); % Number of second planet (arrival)
        
        % Propagation
        [sds,error]=kepPro([xp1,vout],tof*alpha*86400,sys.mu); % Keplerian propagation
        xds=sds(1:3);vds1=sds(4:6);
        if error
            dv=1e6;
            vin=[];
            return
        end
        if options(5) && nargout > 4 % Output requested
            if alpha > 0 % There is a DSM, so this leg exist
                temp = car2kep([xp1,vout], sys.mu);
                output.leg(sindex-1).part(1).kep5 = temp(1:5);
                output.leg(sindex-1).part(1).theta1 = temp(6); % True anomaly at the beginning of the arc
                temp = car2kep(sds, sys.mu);
                output.leg(sindex-1).part(1).theta2 = temp(6); % True anomaly at the end of the arc
                output.leg(sindex-1).part(1).P = 2*pi*sqrt(temp(1)^3/sys.mu) / 86400; % Period [d]
                output.leg(sindex-1).part(1).t1 = time;
                output.leg(sindex-1).part(1).t2 = time+tof*alpha;
                output.leg(sindex-1).part(1).v1 = vout;
                output.leg(sindex-1).part(1).v2 = vds1;
                output.leg(sindex-1).part(1).x1 = xp1;
                output.leg(sindex-1).part(1).x2 = xds;
                output.leg(sindex-1).part(1).v1_rel_xyz = vout - vp1;
                output.leg(sindex-1).part(1).v1_rel_rth = car_rthT(output.leg(sindex-1).part(1).v1_rel_xyz, [xp1,vp1]);
                output.leg(sindex-1).part(1).body1 = planet1;
                output.leg(sindex-1).part(1).v2_rel_xyz = [];
                output.leg(sindex-1).part(1).v2_rel_rth = [];
                output.leg(sindex-1).part(1).body2 = [];
            end
        end
        
        % Lambert arc
        [xp2,vp2] = sys.eph_handle(planet2,time+tof); % Position and velocity of the arrival planet
        [a_lambert,p_lambert,e_lambert,error,vds2,vin] = lambertMR(xds,xp2,tof*(1-alpha)*86400,sys.mu,multirev(1,sindex-1),multirev(2,sindex-1),multirev(3,sindex-1),0); %#ok
        
        if error
            dv=1e6;
            vin=[];
            return
        end
        if options(5) && nargout > 4 % Output requested
            if alpha == 0 % Virtually no DSM, so the Lambert arc is the only physical leg
                part = 1;
                output.leg(sindex-1).part(part).body1 = planet1;
                output.leg(sindex-1).part(part).v1_rel_xyz = vds2 - vp1;
                output.leg(sindex-1).part(part).v1_rel_rth = car_rthT(output.leg(sindex-1).part(part).v1_rel_xyz, [xp1,vp1]);
            else
                part = 2;
                output.leg(sindex-1).part(part).body1 = [];
                output.leg(sindex-1).part(part).v1_rel_xyz = [];
                output.leg(sindex-1).part(part).v1_rel_rth = [];
            end
            temp = car2kep([xds, vds2], sys.mu);
            output.leg(sindex-1).part(part).kep5 = temp(1:5);
            output.leg(sindex-1).part(part).theta1 = temp(6); % True anomaly at the beginning of the arc
            temp = car2kep([xp2, vin], sys.mu);
            output.leg(sindex-1).part(part).theta2 = temp(6); % True anomaly at the end of the arc
            output.leg(sindex-1).part(part).P = 2*pi*sqrt(temp(1)^3/sys.mu) / 86400; % Period [d]
            output.leg(sindex-1).part(part).t1 = time+tof*alpha;
            output.leg(sindex-1).part(part).t2 = time+tof;
            output.leg(sindex-1).part(part).v1 = vds2;
            output.leg(sindex-1).part(part).v2 = vin;
            output.leg(sindex-1).part(part).x1 = xds;
            output.leg(sindex-1).part(part).x2 = xp2;
            output.leg(sindex-1).part(part).v2_rel_xyz = vin - vp2;
            output.leg(sindex-1).part(part).v2_rel_rth = car_rthT(output.leg(sindex-1).part(part).v2_rel_xyz, [xp2,vp2]);
            output.leg(sindex-1).part(part).body2 = planet2;
        end
        
        Tindex=Tindex+2; % Index of vector T
        sindex=sindex+1; % Next planet
        
        if options(4)==2 % Verbose
            fprintf('  ***\nDeep space flight %d-%d. time1 = %f d,MJD2K\n',planet1,planet2,time);
            fprintf('Initial absolute cartesian velocity [km/s]:\n [%f, %f, %f]\n',vout(1),vout(2),vout(3));
            temp = car2rth(vout-vp1, [xp1,vp1]);
            fprintf('Initial relative velocity in r-th-h (also called R-S-T) [km/s]:\n [%f %f %f]\n',temp(1),temp(2),temp(3));
            fprintf('Tof to the next planet = %f d; alpha to DSM = %f\n',tof,alpha);
            temp1 = car2kep([xp1,vout], sys.mu); % Orbital parameters just after the beginning of the leg
            temp2 = car2kep(sds, sys.mu); % Orbital parameters just before the DSM
            fprintf('Orbital parameters before the DSM:\n');
            fprintf(' a = %f km; e = %f; i = %f deg;\n',temp2(1),temp2(2),temp2(3)/pi*180);
            fprintf(' OM = %f deg; om = %f deg; th = %f deg;\n',temp2(4)/pi*180,temp2(5)/pi*180,temp2(6)/pi*180);
            fprintf(' rp = %f km; ra = %f km; P = %f d.\n',temp2(1)*(1-temp2(2)),temp2(1)*(1+temp2(2)),2*pi*sqrt(temp2(1)^3/sys.mu)/86400);
            nr1_time = (tof*alpha) / (2*pi*sqrt(temp1(1)^3/sys.mu)/86400); % Number of revolutions computed in time
            nr1 = floor((tof*alpha) / (2*pi*sqrt(temp1(1)^3/sys.mu)/86400)); % Number of complete revolutions
            if temp2(6)<temp1(6)
                nr1 = nr1 + (2*pi-(temp1(6) - temp2(6)))/(2*pi);
            else
                nr1 = nr1 + (temp2(6) - temp1(6))/(2*pi);
            end
            fprintf('DSM time = %f d,MJD2K; dv = %f km/s\n',tof*alpha,norm(vds2-vds1));
            fprintf('DSM cartesian position [km]:\n [%f, %f, %f]\n',xds(1),xds(2),xds(3));
            temp1 = car2kep([xds, vds2], sys.mu); % Orbital parameters just after the DSM
            temp2 = car2kep([xp2, vin], sys.mu); % Orbital parameters just before the end of the leg
            fprintf('Orbital parameters after the DSM:\n');
            fprintf(' a = %f km; e = %f; i = %f deg;\n',temp1(1),temp1(2),temp1(3)/pi*180);
            fprintf(' OM = %f deg; om = %f deg; th = %f deg;\n',temp1(4)/pi*180,temp1(5)/pi*180,temp1(6)/pi*180);
            fprintf(' rp = %f km; ra = %f km.;\n',temp1(1)*(1-temp1(2)),temp1(1)*(1+temp1(2)));
            fprintf('Final absolute cartesian velocity [km/s]:\n [%f, %f, %f]\n',vin(1),vin(2),vin(3));
            nr2_time = (tof*(1-alpha)) / (2*pi*sqrt(temp1(1)^3/sys.mu)/86400);
            nr2 = floor((tof*(1-alpha)) / (2*pi*sqrt(temp1(1)^3/sys.mu)/86400)); % Number of complete revolutions
            if temp2(6)<temp1(6)
                nr2 = nr2 + (2*pi-(temp1(6) - temp2(6)))/(2*pi);
            else
                nr2 = nr2 + (temp2(6) - temp1(6))/(2*pi);
            end
            temp = car2rth(vin-vp2,[xp2,vp2]);
            fprintf('Final relative velocity in r-th-h (also called R-S-T) [km/s]:\n [%f %f %f]\n',temp(1),temp(2),temp(3));
            temp = vin - vp2;
            fprintf('Final relative velocity in x-y-z [km/s]:\n [%f %f %f]\n', temp(1), temp(2), temp(3));
            fprintf('Lambert parameters: (0=D, 1=R) = %d; nrev = %d; Energy = %d\n',multirev(1,sindex-2),multirev(2,sindex-2),multirev(3,sindex-2));
            fprintf('No. of revolutions on this leg (propagation + lambert):\n %f fractions of P; %f fractions of 2*pi.\n', nr1_time+nr2_time, nr1+nr2);
            if planet1 == planet2 % Resonant leg
                kepp = sys.eph_kep_handle(planet1, time); % Keplerian parameters of the planet
                fprintf('Resonant leg. Resonance: %f:%f (spacecraft:planet).\n', nr1_time+nr2_time, tof/(2*pi*sqrt(kepp(1)^3/sys.mu)/86400));
            end
        end

        % Plot
        if options(4) % Plot flag
            % Plots the spacecraft trajectory from this leg
            if alpha>0
                %[tsp, ysp]=ode45(@two_body_dynamics,[0, tof*alpha*86400],[xp1,vout],ode45options, sys.mu); %#ok
                [tsp, ysp]=ode45(@dyn_2BP, [0, tof*alpha*86400],[xp1,vout],ode45options, sys.mu); %#ok
                line(ysp(:,1), ysp(:,2), ysp(:,3), 'Color', 'k');
            end
            %[tsp, ysp]=ode45(@two_body_dynamics,[0, tof*(1-alpha)*86400],[xds,vds2],ode45options, sys.mu); %#ok
            [tsp, ysp]=ode45(@dyn_2BP,[0, tof*(1-alpha)*86400],[xds,vds2],ode45options, sys.mu); %#ok
            line(ysp(:,1),ysp(:,2),ysp(:,3),'Color','b');
            line(ysp(1,1),ysp(1,2),ysp(1,3),'Marker','.') % Marker in dsm point
            line(xp2(1),xp2(2),xp2(3),'Marker','o') % Marker in rendevouz with planet point
        end
        
        % Calculation of dv of deep space manoeuvre
        dv = [dv;norm(vds2-vds1)];
        % Initialises variables for next flyby
        time=time+tof; % Time of the next flyby (during istantaneous flyby)
        % vin (velocity vector at the next planet encounter) already initialised
        planet1=planet2; % Next planet
        xp1=xp2; % Position of the next planet
        vp1=vp2; % Velocity of the next planet
    end
end

if options(4) % Plot flag
    % Plots all the planets from t0 to (t0 + sum(tof))=time
    s_un = unique(s);
    a_color = s_un; % Vector of colours for each planet in the sequence
    % Normalisation: a_color between 0 and 1
    a_color = a_color - min(a_color);
    if ~all(a_color == 0)
        a_color = a_color ./ max(a_color);
    end
    for i=1:length(s_un)
        kepp = sys.eph_kep_handle(s_un(i), time); % Keplerian parameters of the planet
        Pp = 2*pi*sqrt(kepp(1)^3/sys.mu) / 86400; % Period of the planet [d]
        step = Pp/360; % Roughly one point for each degree
        ti = [t0:step:time, time]; % Vector of times [d]
        xpi = zeros(length(ti), 3);
        for k = 1:length(ti)
            xpi(k,:) = sys.eph_handle(s_un(i), ti(k)); % Position and velocity of the arrival planet % Position and velocity of the arrival planet
        end
        line(xpi(:,1), xpi(:,2), xpi(:,3), 'Color',[mod(a_color(i),1.001), mod(a_color(i)+.33,1.001), mod(a_color(i)+0.66,1.001)]);
    end
    
    % Finalising plot
    axis equal;
    hold off;
end