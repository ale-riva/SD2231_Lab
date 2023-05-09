function x_n = Vehicle_state_eq(x,delta)
% ADDME Dynamic model function
%    x = the states from the sigma points (yawrate in rad)
%    param = parameters that you might need, such as vehicle parameters.
%    delta = steering angle (in rad) for previous time instance

global lf lr mass Iz Cf Cr
vx = x(1);
vy = x(2);
yawrate = x(3);
dt = 0.01;

alpha12 = atan((vy+yawrate*lf)/vx)-delta;
alpha34 = atan((vy-yawrate*lr)/vx);
F12 = -Cf*alpha12;
F34 = -Cr*alpha34;
  

vx_dot = -(F12*sin(delta))/(mass)+yawrate*vy;
vy_dot = ((F34+F12*cos(delta))/(mass))-yawrate*vx;
yawrate_dot = (lf*F12*cos(delta)-lr*F34)/(Iz);

f_x = [vx_dot;vy_dot;yawrate_dot];


% Integrate using Runge Kutta (in the script folder) or simple euler forward

f = @(x)[f_x(1,:);f_x(2,:);f_x(3,:)];
x_n = rk4(f,dt,x(1:3,:));