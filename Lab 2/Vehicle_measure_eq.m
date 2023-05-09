function y_n = Vehicle_measure_eq(x,delta)
% ADDME Measurement function
%    x = the states
%    param = parameters that you might need, such as vehicle parameters.

global lf lr mass Cf Cr

vx = x(1);
vy = x(2);
yawrate = x(3);

alpha12 = atan2((vy+yawrate*lf),vx)-delta;
alpha34 = atan2((vy-yawrate*lr),vx);
F12 = -Cf*alpha12;
F34 = -Cr*alpha34;

ay = (F34+F12*cos(delta))/(mass);
y_n = [vx;ay;yawrate];
