function y_n = Vehicle_measure_eq(x,delta)
% ADDME Measurement function
%    x = the states
%    param = parameters that you might need, such as vehicle parameters.

global lf lr mass Mu Cf Cr Ly_relax

vx = x(1);
vy = x(2);
yawrate = x(3);

alpha12 = atan((vy+yawrate*lf)/vx)-delta;
alpha34 = atan((vy-yawrate*lr)/vx);

lambda12 = (Mu*0.5*mass*9.82)/(2*Cf*abs(tan(alpha12)));
lambda34 = (Mu*0.5*mass*9.82)/(2*Cr*abs(tan(alpha34)));

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

ay = (F34_dug+F12_dug*cos(delta))/(mass);
y_n = [vx;ay;yawrate];
