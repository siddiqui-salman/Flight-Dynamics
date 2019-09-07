%Team Gemini
%Project Final Deliverable

clear all
close all
clc


%%   Flight Condition
inialt = 55000;                      %   altitude at trim (ft)
rho = .002378*(1-.0000068756*...%   air density at altitude (kg/ft^3)
    inialt)^4.2561;                 
T = 518.67*(1-.0000068756*inialt);   %   Temperature at altitude (deg R)    
a = 967.78           %   speed of sound at altitude (ft/s)

%%   Trim Conditions
%   Trim Linear Velocity 

Vt0 = 1742;  %   trim total velocity (ft/s)
M_T = 1.8;
alpha0 = 2.9737954*pi/180;             %   trim angle of attack (rad)
beta0 = 0;                          %   trim sidelip angle (rad)
u0 = Vt0*cos(alpha0)*cos(beta0);    %   trim body-x (forward) velocity (ft/s)
v0 = Vt0*sin(beta0);                %   trim body-y (side) velocity (ft/s)
w0 = Vt0*sin(alpha0)*cos(beta0);    %   trim body-z velocity (ft/s)

Xw = 2*[52*Vt0, 104*Vt0, 104*Vt0, 52*Vt0, 0];
Yw = [30*Vt0, -30*Vt0, 30*Vt0, -30*Vt0, 0];
%   Trim Angular Velocity (default = trim values)

p0 = 0;                      %   trim body-x (roll) angular velocity (rad/s)
q0 = 0;                      %   trim body-y (pitch) angular velocity (rad/s)
r0 = 0;                      %   trim body-z (yaw) angular velocity (rad/s)

%   Trim Attitude (roll, pitch)

phi0 = 0;                   %   trim roll angle (rad)
theta0 = alpha0;            %   trim pitch angle (rad)

%   Trim Control Inputs 

de0 = -4.1938979*pi/180;            %   trim elevator deflection (rad)
dT0 = 0.4868;                   %   trim throttle setting (% of full power)
da0 = 0;                        %   trim aileron deflection (rad)
dr0 = 0;                        %   trim rudder deflection (rad)

%%   Initial Conditions
%   Initial Linear Velocity (default values are trim values)

initu = u0;                     %   initial body-x (forward) velocity (ft/s)
initv = v0;                     %   initial body-y (side) velocity (ft/s)
initw = w0;                     %   initial body-z velocity (ft/s)

%   Initial Angular Velocity (default = trim values)

initp = p0;                     %   initial body-x (roll) angular velocity (rad/s)
initq = q0;                     %   initial body-y (pitch) angular velocity (rad/s)
initr = r0;                     %   initial body-z (yaw) angular velocity (rad/s)

%   Initial Attitude (roll, pitch, yaw)

initroll = phi0;                %   initial roll angle (rad)
initpitch = theta0;             %   initial pitch angle (rad)
inityaw = 0.0;                  %   initial yaw angle (rad)
    
%   Initial Position 

initnorth = 0;                  %   initial north position (ft)
initeast = 0;                   %   initial east position (ft)
initalt = inialt;                    %   initial altitude (ft)

%%   Starfighter aerodynamic coefficients

%   CL (Lift Coefficient)

cla = 2.005;               %   change in CL with respect to alpha (angle of attack) (/rad)
cladot = 0.82;             %   change in CL with respect to normalized alphadot (/rad)
clq = 1.9;              %   change in CL with respect to normalized pitch rate (/rad)
clde = 0.523;              %   change in CL with respect to elevator deflection (/rad)
clo = 0.122;               %   CL at zero angle of attack 

%   CD (Drag Coefficient) 

cda = 0.384;             %   change in CD with respect to alpha (angle of attack) (/rad)
cdo = 0.048;             %   CD at zero angle of attack 
cdde = 0;               %   change in CD with respect to elevator deflection (/rad)

%   CM (Pitch Moment Coefficient) 

cmo = -0.028;            %   CM at zero angle of attack
cma = -1.308;             %   change in CM with respect to alpha (angle of attack) (/rad)
cmadot = -2.05;         %   change in CM with respect to normalized alphadot (/rad)
cmq = -4.83;            %   change in CM with respect to normalized pitch rate (/rad)
cmde = -1.31;           %   change in CM with respect to elevator deflection (/rad)

%   CY (Side Force Coefficient)

cyb = -1.045;           %   change in CY with respect to sideslip (/rad)
cyp = 0;          %   change in CY with respect to normalized roll rate (/rad)
cyr = 0;              %   change in CY with respect to normalized yaw rate (/rad)
cyda = 0;               %   change in CY with respect to aileron deflection (/rad)
cydr = 0.087;             %   change in CY with respect to rudder deflection (/rad)    

%   CL (Roll Moment Coefficient)

clb = -0.093;          %   change in Cl with respect to sideslip (/rad)
clp = -0.272;           %   change in Cl with respect to normalized roll rate (/rad)
clr = 0.154;           %   change in Cl with respect to normalized yaw rate (/rad)
clda = 0.0173;           %   change in Cl with respect to aileron deflection (/rad)
cldr = 0.0079;           %   change in Cl with respect to rudder deflection (/rad)

%   CN (Yaw Moment Coefficient)

cnb = 0.242;           %   change in CN with respect to sideslip (/rad)
cnp = -0.093;          %   change in CN with respect to normalized roll rate (/rad)
cnr = -0.649;           %   change in CN with respect to normalized yaw rate (/rad)
cnda = -0.0025;         %   change in CN with respect to aileron deflection (/rad)
cndr = -0.0435;         %   change in CN with respect to rudder deflection (/rad)

%%   Starfighter Mass and Geometry Properties

sw = 196;               %   wing surface area (ft^2)
b = 21.9;               %   wing span (ft)
cbar = 9.6;            %   mean chord of wing (ft)
weight = 16300;          %   weight (lbs)

T_max = 11905;           %   maximum engine thrust (lbs)
cg = [0 0 0];           %   location of CG in body coordinates (ft)
ac = [0 0 0];           %   location of AC in body coordinates (ft)
eng = [0 0 0];          %   location of engine in body coordinates (ft)

Ix = 3600;              %   moment of inertia about the body-x axis (slug*ft^2)
Iy = 59000;              %   moment of inertia about the body-y axis (slug*ft^2)
Iz = 60000;              %   moment of inertia about the body-z axis (slug*ft^2)
Ixz = 0;                %   product of inertia about the body-xz axes (slug*ft^2)

%%   Waypoint Parameters
K_phi =10; % Update

%%   Run Simulation
dt = 0.01;
time = 0:dt:2500;

[time,x,y] = sim('Starfighter_6DoF_runsim_waypoint',max(time)); 

%   Plot Results

u = y(:,1);                 %   Body-X Velocity (ft/s)
v = y(:,2);                 %   Body-Y Velocity (ft/s)
w = y(:,3);                 %   Body-Z Velocity (ft/s)

p = 180/pi*y(:,4);          %   Body-X Angular Velocity (deg/s)
q = 180/pi*y(:,5);          %   Body-Y Angular Velocity (deg/s)
r = 180/pi*y(:,6);          %   Body-Z Angular Velocity (deg/s)

phi = 180/pi*y(:,7);        %   Roll Angle (deg)
theta = 180/pi*y(:,8);      %   Pitch Angle (deg)
psi = 180/pi*y(:,9);        %   Yaw Angle (deg)

X_N = y(:,10);              %   North Position (ft)
Y_E = y(:,11);              %   East Position (ft)
inialtp = y(:,12);                %   Altitude (ft)

d_e = 180/pi*y(:,13);       %   Elevator Deflection (deg)
d_r = 180/pi*y(:,15);       %   Rudder Deflection (deg)
d_a = 180/pi*y(:,14);       %   Aileron Deflection (deg)
d_T = y(:,16);              %   Throttle Control

%   Translational Velocity

figure
subplot(311)
plot(time,u)
title('Translational Velocity')
ylabel('u (ft/s)')
subplot(312)
plot(time,v)
ylabel('v (ft/s)')
subplot(313)
plot(time,w)
ylabel('w (ft/s)')
xlabel('Time (s)')

%   Angular Velocity

figure(2)
subplot(311)
plot(time,p)
title('Angular Velocity')
ylabel('p (deg/s)')
subplot(312)
plot(time,q)
ylabel('q (deg/s)')
subplot(313)
plot(time,r)
ylabel('r (deg/s)')
xlabel('Time (s)')

%   Euler Angles (Orientation)

figure(3)
subplot(311)
plot(time,phi)
title('Attitude (Euler Angles)')
ylabel('Roll (deg)')
subplot(312)
plot(time,theta)
ylabel('Pitch (deg)')
subplot(313)
plot(time,psi)
ylabel('Yaw (deg)')
xlabel('Time (s)')

%   Inertial Position

figure(4)
subplot(311)
plot(time,X_N)
title('Inertial (NED) Position')
ylabel('North (ft)')
subplot(312)
plot(time,Y_E)
ylabel('East (ft)')
subplot(313)
plot(time,inialtp)
ylabel('Altitude (ft)')
xlabel('Time (s)')

%   Control Inputs

figure(5)
subplot(411)
plot(time,d_e)
title('Control Inputs')
ylabel('Elevator (deg)')
subplot(412)
plot(time,d_T)
ylabel('Throttle')
subplot(413)
plot(time,d_a)
ylabel('Aileron (deg)')
subplot(414)
plot(time,d_a)
ylabel('Rudder (deg)')
xlabel('Time (s)')

figure
hold
title('Horizontal Flight Path')
plot(Y_E,X_N,Yw,Xw,'r*')
xlabel('East (ft)')
ylabel('North (ft)')
hold

figure
subplot(211)
plot(time,y(:,17))
title('psi_b')
subplot(212)
plot(time,y(:,18))
title('psi')