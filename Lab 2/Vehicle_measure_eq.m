function y_n = Vehicle_measure_eq(x,delta)
% ADDME Measurement function
%    x = the states
%    param = parameters that you might need, such as vehicle parameters.

global lf lr mass Cf Cr

vx = x(1);
vy = x(2);
yawrate = x(3);
% disp(size(x))

F12 = -Cf*atan((vy+yawrate*lf)/vx)-delta;
F34 = -Cr*atan((vy-yawrate*lr)/vx);

ay = (F34+F12*cos(delta))/(mass);
y_n = [vx;ay;yawrate];
