%Team Gemini
%Project Final Deliverable
%Root Locus plots.


clc;
clear;

h = 55000;                          %   altitude at trim (ft)
rho = .002378*(1-.0000068756*...    %   air density at altitude (slugs/ft^3)
    h)^4.2561;                 
T = 518.67*(1-.0000068756*h);       %   Temperature at altitude (deg R)    
a = 967.78;               %   speed of sound at altitude (ft/s)
g = 32.174;                         %   gravity
Vt0 = 1742;                      %   trim total velocity (ft/s)
Mach0 = 1.8;                      %   trim mach number
alpha0 = 2.9737954*pi/180;             %   trim angle of attack (rad)
de0 =  -4.1938979*pi/180;                        %   elevator trim
dT0 =  5795.5876705721688901230688059874;                     %   thrust trim
beta0 = 0;                          %   trim sidelip angle (rad)
u0 = Vt0*cos(alpha0)*cos(beta0);    %   trim body-x (forward) velocity (ft/s)
v0 = Vt0*sin(beta0);                %   trim body-y (side) velocity (ft/s)
w0 = Vt0*sin(alpha0)*cos(beta0);    %   trim body-z velocity (ft/s)

%   CL (Lift Coefficient)

cla = 2.005;               %   change in CL with respect to alpha (angle of attack) (/rad)
cladot = 0.82;             %   change in CL with respect to normalized alphadot (/rad)
clq = 1.9;              %   change in CL with respect to normalized pitch rate (/rad)
clde = 0.523;              %   change in CL with respect to elevator deflection (/rad)
clo = 0.122;               %   CL at zero angle of attack 
cl_mach =-0.111;              %   change in CL with respect to mach number
cl_trim = 0.18778;         %   cl trim

%   CD (Drag Coefficient) 

cda = 0.384;             %   change in CD with respect to alpha (angle of attack) (/rad)
cdo = 0.048;             %   CD at zero angle of attack 
cdde = 0;               %   change in CD with respect to elevator deflection (/rad)
cd_mach = 0;              %   change in Cd with respect to mach number
cd_trim = 0.432;           % cd trim

%   CM (Pitch Moment Coefficient) 

cmo = -0.028;            %   CM at zero angle of attack
cma =  -1.308;             %   change in CM with respect to alpha (angle of attack) (/rad)
cmadot = -2.05;         %   change in CM with respect to normalized alphadot (/rad)
cmq = -4.83;            %   change in CM with respect to normalized pitch rate (/rad)
cmde =  -1.31;           %   change in CM with respect to elevator deflection (/rad)
cm_mach = -0.033;            %   change in Cm with respect to mach number
cm_trim =0;            % cm trim

%   CY (Side Force Coefficient)

cyb = -1.045;           %   change in CY with respect to sideslip (/rad)
cyp =  0;          %   change in CY with respect to normalized roll rate (/rad)
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

cnb =0.242;           %   change in CN with respect to sideslip (/rad)
cnp =-0.093;          %   change in CN with respect to normalized roll rate (/rad)
cnr = -0.649;           %   change in CN with respect to normalized yaw rate (/rad)
cnda =  -0.0025;         %   change in CN with respect to aileron deflection (/rad)
cndr =  -0.0435;         %   change in CN with respect to rudder deflection (/rad)

%%   Starfighter Mass and Geometry Properties

sw = 196;               %   wing surface area (ft^2)
b = 21.9;               %   wing span (ft)
cbar = 9.6;            %   mean chord of wing (ft)
m = 16300/g;          %   weight (lbs)
qw = 0.5*rho*Vt0^2;

Ix = 3600;              %   moment of inertia about the body-x axis (slug*ft^2)
Iy =  59000;              %   moment of inertia about the body-y axis (slug*ft^2)
Iz = 60000;              %   moment of inertia about the body-z axis (slug*ft^2)
Ixz = 0;                %   product of inertia about the body-xz axes (slug*ft^2)
cg = [0 0 0];           %   location of CG in body coordinates (ft)
ac = [0 0 0];           %   location of AC in body coordinates (ft)
eng = [0 0 0];          %   location of engine in body coordinates (ft)
%   Trim Angular Velocity (default = trim values)

p0 = 0;                      %   trim body-x (roll) angular velocity (rad/s)
q0 = 0;                      %   trim body-y (pitch) angular velocity (rad/s)
r0 = 0;                      %   trim body-z (yaw) angular velocity (rad/s)

%   Trim Attitude (roll, pitch)

phi0 = 0;                   %   trim roll angle (rad)
theta0 = alpha0;            %   trim pitch angle (rad)

%   Trim Control Inputs 


%dT0 = 0.48;                   %   trim throttle setting (% of full power)
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
initalt = 0;                    %   initial altitude (ft)


%% Longitudnal stability derivativs
Du = (qw*sw/u0)*(2*cd_trim+Mach0*cd_mach);
Da = (qw*sw)*(cda-cl_trim);
Dadot =0;
Dq = 0;
Dde = qw*sw*cdde;
Lu = (qw*sw/u0)*(2*cl_trim+Mach0*cl_mach);
La = (qw*sw)*(cla+cd_trim);
Ladot = ((qw*cbar*sw)/(u0*2))*cladot;
Lq = ((qw*cbar*sw)/(u0*2))*clq;
Lde = qw*sw*clde;
Mu = ((qw*cbar*sw*Mach0)/(u0))*cm_mach;
Ma = (qw*cbar*sw)*cma;
Mw = (1/u0)*Ma;
Madot = ((qw*cbar^2*sw)/(u0*2))*cmadot; 
Mwdot = (1/u0)*Madot;
Mq = ((qw*cbar^2*sw)/(u0*2))*cmq;
Mde = (qw*cbar*sw)*cmde;

%% Lateral stability derivitives
Yb = (qw*sw)*cyb;
Yp =((qw*b*sw)/(u0*2))*cyp;
Yv = (1/u0)*Yb;
Yr = ((qw*b*sw)/(u0*2))*cyr;
Yda = (qw*sw)*cyda;
Ydr = (qw*sw)*cydr;
Lb = (qw*sw*b)*clb;
Lp =((qw*b^2*sw)/(u0*2))*clp;
Lv = (1/u0)*Lb;
Lr = ((qw*b^2*sw)/(u0*2))*clr;
Lda = (qw*sw*b)*clda;
Ldr = (qw*sw*b)*cldr;
Nb = (qw*sw*b)*cnb;
Np =((qw*b^2*sw)/(u0*2))*cnp;
Nv = (1/u0)*Nb;
Nr = ((qw*b^2*sw)/(u0*2))*cnr;
Nda = (qw*sw*b)*cnda;
Ndr = (qw*sw*b)*cndr;

%% longitudnal derivatives
Xu = -Du*cos(alpha0) + Lu*sin(alpha0);
Xa = -Da*cos(alpha0) + La*sin(alpha0);
Xw = (1/u0)*Xa;
Xq = -Dq*cos(alpha0) + Lq*sin(alpha0);
Xde = -Dde*cos(alpha0) + Lde*sin(alpha0);
Xadot = -Dadot*cos(alpha0) + Ladot*sin(alpha0);
Xwdot = (1/u0)*(Xadot);

Zu = -Du*sin(alpha0) - Lu*cos(alpha0) ;
Za = -Da*sin(alpha0) - La*cos(alpha0) ;
Zw = (1/u0)*Za;
Zq = -Dq*sin(alpha0) - Lq*cos(alpha0) ;
Zde = -Dde*sin(alpha0) - Lde*cos(alpha0) ;
Zadot = -Dadot*sin(alpha0) - Ladot*cos(alpha0) ;
Zwdot = (1/u0)*Zadot;

TdT = 11905;      % Tmax = dT/dT
%% Longitudnal Matrix
a11 = (Xu/m);
a12 = (Xw/m);
a13 = (Xq/m-w0);
a14 = -g*cos(alpha0);
a21 = Zu/(m-Zwdot);
a22 = Zw/(m-Zwdot);
a23 = (Zq+m*u0)/(m-Zwdot);
a24 = (-m*g*sin(alpha0))/(m-Zwdot);
a31 = (1/Iy)*(Mu+((Zu*Mwdot)/(m-Zwdot)));
a32 = (1/Iy)*(Mw+((Zw*Mwdot)/(m-Zwdot)));
a33 = (1/Iy)*(Mq+(((Zq+m*u0)*Mwdot)/(m-Zwdot)));
a34 = (-1/Iy)*(Mwdot*m*g*sin(alpha0))*(1/(m-Zwdot));
a41 = 0;
a42 = 0;
a43 =1;
a44 =0;

b11 = (Xde/m);
b12 = (TdT/m);
b21 = Zde/(m-Zwdot);
b22 = 0;
b31 = (1/Iy)*(Mde+((Zde*Mwdot)/(m-Zwdot)));
b32 = 0;
b41 = 0;
b42 = 0;
theta0 = alpha0;

A_long = [a11 a12 a13 a14 0;a21 a22 a23 a24 0;a31 a32 a33 a34 0;a41 a42 a43 a44 0;sin(theta0) -cos(theta0) 0 Vt0 0];
B_long = [b11 b12;b21 b22;b31 b32;b41 b42;0 0];


[eigvec_long,Deigval_long] = eig(A_long);
lam_long = diag(Deigval_long);

%% Lateral Matrix
I1 = Ix*Iz-Ixz^2;

c11 =Yv/m;
c12 =(Yp/m+w0);
c13 =(Yr/m-u0);
c14 = g*cos(alpha0);
c21 = (Iz/I1)*Lv+(Ixz/I1)*Nv;
c22 = (Iz/I1)*Lp+(Ixz/I1)*Np;
c23 = (Iz/I1)*Lr+(Ixz/I1)*Nr;
c24 =0;
c31 =(Ix/I1)*Nv+(Ixz/I1)*Lv;
c32 =(Ix/I1)*Np+(Ixz/I1)*Lp;
c33 =(Ix/I1)*Nr+(Ixz/I1)*Lr;
c34 =0;
c41 =0;
c42 =1;
c43 =tan(alpha0);
c44 =0;

d11 = (Yda/m);
d12 = (Ydr/m);
d21 = (Iz/I1)*Lda+(Ixz/I1)*Nda;
d22 = (Iz/I1)*Ldr+(Ixz/I1)*Ndr;
d31 = (Ix/I1)*Nda+(Ixz/I1)*Lda;
d32 = (Ix/I1)*Ndr+(Ixz/I1)*Ldr;
d41 = 0;
d42 = 0;

A_lat = [c11 c12 c13 c14 0;c21 c22 c23 c24 0; c31 c32 c33 c34 0; c41 c42 c43 c44 0;0 0 sec(theta0) 0 0];
B_lat = [d11 d12;d21 d22;d31 d32;d41 d42;0 0];


[eigvec_lat,Deigval_lat] = eig(A_lat);
lam_lat = diag(Deigval_lat);

%% Find Transfer Functions

sys_long = ss(A_long,B_long,eye(5),0);
sys_lat = ss(A_lat,B_lat,eye(5),0);


tf_long = tf(sys_long); %[du/de du/dt; dw/de dw/dt; dq/de dq/dt; dtheta/de dtheta/dt;dh/de dh/dt]
tf_lat = tf(sys_lat); %[dv/da dv/dr; dp/da dp/dr; dr/da dr/dr; dphi/da dphi/dr;dpsi/da dpsi/dr]
%% Dimentional Coefficients
% Lateral 
D.Yb = Yb/m;
D.Yr = Yr/m;
D.Ydr = Ydr/m;
D.Nb = Nb/Iz;
D.Nr = Nr/Iz;
D.Ndr = Ndr/Iz;
D.Lp = Lp/Ix;
D.Lda = Lda/Ix;
% Longitudnal
D.Mde = Mde/Iy;
D.Mq = Mq/Iy;
D.Ma= Ma/Iy;
D.Madot =Madot/Iy; 
D.Zde = Zde/m;
D.Za = Za/m;
D.XdT = TdT/m;
D.Zu =Zu/m;

I_S = tf([1],[1 0]);
dtheta_de = tf_long(4,1);
dq_de = tf_long(3,1);
dphi_da = tf_lat(4,1);%%tf([D.Lda],[1 -D.Lp 0]);%%
dbeta_dr = tf_lat(1,2)/Vt0;
dh_de = tf_long(5,1);
du_dt = tf_long(1,2);
dpsi_dr = tf_lat(5,2);

%% Transfer function for Actuator
tf_ele = tf(-10,[1 10]);
tf_ail = tf(10,[1 10]);
tf_rud = tf(10,[1 10]);
tf_rudhead = tf(10,[1 10]);
tf_servo = tf(0.1,[1 0.1]);
tf_actuator = tf(10,[1 10]);
%% TRANSFER FUNCTIONS for Ziegler-Nicholas Method
C_Rcont = tf_ail*dphi_da;
C_betacont = tf_rud*dbeta_dr;
C_altcont = tf_ele*dh_de;
C_velcont = tf_actuator*tf_servo*du_dt;
C_heading = tf_rudhead*dpsi_dr;
%% Apply Ziegler-Nicholas Method

%Roll control
rlocusplot(C_Rcont)
grid on;
pause(20)
[K_roll,poles_roll] = rlocfind(C_Rcont);
Tu_roll = (2*pi)/imag(poles_roll(2));
Kp_roll = 0.6*K_roll;
Ki_roll = (0.6/0.5)*(K_roll/Tu_roll);
Kd_roll = (0.6/0.125)*(K_roll/Tu_roll);

%SIDE SLIP CONTROL
rlocus(C_betacont);
grid on;
pause(10)
[K_beta,poles_beta] = rlocfind(C_betacont);
Tu_betacont = (2*pi)/imag(poles_beta(2));
Kp_beta = 0.6*K_beta;
Ki_beta = (0.6/0.5)*(K_beta/Tu_betacont);
Kd_beta = (0.6/0.125)*(K_beta/Tu_betacont);
% 
% ALTITUDE CONTROL
rlocus(C_altcont);
grid on;
pause(10)
[K_ALT,poles_ALT] = rlocfind(C_altcont);
Tu_ALT = (2*pi)/imag(poles_ALT(2));
Kp_ALT = 0.6*K_ALT;
Ki_ALT = (0.6/0.5)*(K_ALT/Tu_ALT);
Kd_ALT = (0.6/0.125)*(K_ALT/Tu_ALT);

% VELOCITY CONTROL
rlocus(C_velcont);
 grid on;
pause(10)
[K_vel,poles_vel] = rlocfind(C_velcont);
Tu_vel = (2*pi)/imag(poles_vel(4));
Kp_vel = 0.6*K_vel;
Ki_vel = (0.6/0.5)*(K_vel/Tu_vel);
Kd_vel = (0.6/0.125)*(K_vel/Tu_vel);

% Heading CONTROL
rlocus(C_heading);
 grid on;
pause(10)
[K_head,poles_head] = rlocfind(C_heading);
Tu_head = (2*pi)/imag(poles_head(4));
Kp_head = 0.6*K_head;
Ki_head = (0.6/0.5)*(K_head/Tu_head);
Kd_head = (0.6/0.125)*(K_head/Tu_head);


display(C_Rcont)
