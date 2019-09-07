%Team Gemini
%Project Final Deliverable
%Testing of the controls.

clc
clear

phiCom=10*pi/180;

betaCom=1*pi/180;
Vt0 = 1742;
vCom=Vt0*sin(betaCom);

hCom=10;

uCom=10;

dt = 0.01;              %   simulation time step (sec)
tf = 50;              %   final simulation time (sec)
time = 0:dt:tf;         %   simulation time vector (sec)

[time,x,y] = sim('Controls',max(time));

phi=180/pi*y(:,1);
da=180/pi*y(:,2);
B=180/pi*asin(y(:,3)/Vt0);
dr=180/pi*y(:,4);
h=y(:,5);
de=180/pi*y(:,6);
u=y(:,7);
dT=y(:,8);

figure(1)
title('Roll')
subplot(211)
plot(time,phi)
ylabel('phi (deg)')
subplot(212)
plot(time,da)
ylabel('da (deg)')
xlabel('Time (s)')

figure(2)
title('Sideslip')
subplot(211)
plot(time,B)
ylabel('B (deg)')
subplot(212)
plot(time,dr)
ylabel('dr (deg)')
xlabel('Time (s)')

figure(3)
title('Altitude')
subplot(211)
plot(time,h)
ylabel('h (deg)')
subplot(212)
plot(time,de)
ylabel('de (deg)')
xlabel('Time (s)')

figure(4)
title('Forward Velocity')
subplot(211)
plot(time,u)
ylabel('u (deg)')
subplot(212)
plot(time,dT)
ylabel('dT')
xlabel('Time (s)')