function y_n = Vehicle_measure_eq(x,delta)
% ADDME Measurement function
%    x = the states
%    param = parameters that you might need, such as vehicle parameters.

global lf lr mass Cf Cr g Mu L

vx = x(1);
vy = x(2);
yawrate = x(3);

alpha12 = atan2((vy+yawrate*lf),vx)-delta;
alpha34 = atan2((vy-yawrate*lr),vx);
Fz = mass*g;
lambda12 = (Mu*mass*g*lf)/(2*L*Cf*abs(tan(alpha12)));
lambda34 = (Mu*mass*g*lr)/(2*Cr*abs(tan(alpha34)));

if lambda12 <=1
    f_lambda12 = lambda12*(2-lambda12);
else 
    f_lambda12 = 1;
end

if lambda34 <=1
    f_lambda34 = lambda34*(2-lambda34);
else 
    f_lambda34 = 1;
end

F12 = -Cf*tan(alpha12)*f_lambda12;
F34 = -Cr*tan(alpha34)*f_lambda34;
% F12 = -Cf*alpha12;
% F34 = -Cr*alpha34;
ay = (F34+F12*cos(delta))/(mass);
y_n = [vx;ay;yawrate];
