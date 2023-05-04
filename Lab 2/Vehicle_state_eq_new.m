function x_n = Vehicle_state_eq(x,delta)
% ADDME Dynamic model function
%    x = the states from the sigma points (yawrate in rad)
%    param = parameters that you might need, such as vehicle parameters.
%    delta = steering angle (in rad) for previous time instance

global lf lr mass Iz Cf Cr Ly_relax Mu
vx = x(1);
vy = x(2);
yawrate = x(3);
% F12 = x(4);
% F34 = x(5);
dt = 0.01;

alpha12 = atan((vy+yawrate*lf)/vx)-delta;
alpha34 = atan((vy-yawrate*lr)/vx);

num = Mu*0.5*mass*9.82;
frac12 = 2*Cf*abs(tan(alpha12));
frac34 = 2*Cf*abs(tan(alpha34));
lambda12 = num/frac12;
lambda34 = num/frac34;

if lambda12 < 1
    slip_func12 = (2-lambda12)*lambda12;
else
    slip_func12 = 1;
end
if lambda34 < 1
    slip_func34 = (2-lambda34)*lambda34;
else
    slip_func34 = 1;
end


F12_dug = -Cf*tan(alpha12)*slip_func12;
F34_dug = -Cr*tan(alpha34)*slip_func34;
  

vx_dot = -(F12_dug*sin(delta))/(mass)+yawrate*vy;
vy_dot = ((F34_dug+F12_dug*cos(delta))/(mass))-yawrate*vx;
yawrate_dot = (lf*F12_dug*cos(delta)-lr*F34_dug)/(Iz);

f_x = [vx_dot;vy_dot;yawrate_dot];

% Integrate using Runge Kutta (in the script folder) or simple euler forward

f = @(x)[f_x(1,:);f_x(2,:);f_x(3,:)];
x_n = rk4(f,dt,x(1:3,:));