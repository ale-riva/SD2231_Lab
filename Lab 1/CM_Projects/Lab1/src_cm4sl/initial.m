    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  KTH SD2231 - Applied Vehicle Dynamics Control
%
%  Lab:     1 - Slip Control for vehicles
%  Date:    Spring term 2020
%  Teacher: Mikael Nybacka and Wenliang Zhang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
clc
close all
clear all

global Veh
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Choose input parameters
mu_select = 1;              % Here you can change to different friction levels
                            % but for this lab you should stick to mu_select = 1
                            % representing dry road. 
                            % If you have time you can play around with
                            % other mu levels but the tyre model in
                            % CarMaker will not be realistic for snow or
                            % wet road.
dt = 0.001;
T = 0.039; L  = 0.052;
T_new_rl_b = 0.287; L_new_rl_b = 0.048;
T_new = 0.0155; L_new = 0.0475;



%final values morning 13/04
% Kp_br_r = 2*T/L;
% Ki_br_r = 2*T/(L^2);
% Kd_br_r = 2*T;
%final values afternoon 13/04
% Kp_br_r = 2*T/L
% Ki_br_r = 0.7*T/(L^2);
% Kd_br_r = 2*T

% Kp_br_r = 1.2*T_new/L_new;
% Ki_br_r= 0.6*T_new/L_new^2;
% Kd_br_r = 0.6*T_new;


%Ziegler-Nichols values 
Kp_tr = 2.2*1.2*T_new/L_new;
Ki_tr = 0.6*T_new/L_new^2;
Kd_tr = 0.6*T_new;

Kp_br_r = 1.2*T_new_rl_b/L_new_rl_b;
Ki_br_r = 0.6*T_new_rl_b/L_new_rl_b^2;
Kd_br_r = 0.6*T_new_rl_b;

Kp_br_f = 1.2*T_new_rl_b/L_new_rl_b;
Ki_br_f = 0.5*0.6*T_new_rl_b/L_new_rl_b^2;
Kd_br_f = 0*0.6*T_new_rl_b;





% Kp_tr =1*T/L;
% Ki_tr =0.4*T/(L^2);
% Kd_tr = 2*T;

N = 100;
delay = 250;


K_awup = Ki_br_r/Kp_br_r;
ref_tr = 0.06;
ref_br_f = 0.15;
ref_br_r = 0.17-0.1; 
threshold = 0.01;
threshold_tc = ref_tr-0.02;
threshold_br_f = ref_br_f-0.05;
threshold_br_r = ref_br_r-0.05;
step_time = 2;
vx_threshold = 1;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vehicle parameters DO NOT CHANGE
% ______________________________________________________________________________
% 
% BASIC PARAMETERS         UNIT          TOTAL
% ______________________________________________________________________________
% 
% Total mass               kg           2107.996
% Total weight             N           20671.009
% Total roll inertia       kgm^2         952.228
% Total pitch inertia      kgm^2        3508.505
% Total yaw inertia        kgm^2        3954.288
% Total CoG height (FrD)   m               0.545
% 
% Sprung mass              kg           1841.860
% Sprung weight            N           18061.279
% Sprung roll inertia      kgm^2         741.760
% Sprung pitch inertia     kgm^2        2903.000
% Sprung yaw inertia       kgm^2        3170.600
% Sprung CoG height (FrD)  m               0.578
% 
% Unsprung mass (total)    kg            266.136
% 
% Wheelbase                m               2.971
% 
% ______________________________________________________________________________
% 
% BASIC PARAMETERS         UNIT      AVERAGE/TOTAL       LEFT         RIGHT
% ______________________________________________________________________________
% 
% Front axle:
% 
% Track width              m               1.662
% Toe Angle                deg             0.200          0.200         0.200
% Camber Angle             deg            -0.000         -0.000        -0.000
% 
% Caster Angle             deg             1.298          1.298         1.298
% Caster Offset            mm              6.077          6.077         6.077
% Kingpin Angle            deg            13.257         13.257        13.257
% Kingpin Offset           mm              2.809          2.809         2.809
% 
% Normal force             N           10436.920       5218.459      5218.461
% Loaded tire radius       m               0.296          0.296         0.296
% Tire rate                N/mm          236.700        236.700       236.700
% ______________________________________________________________________________
% 
% Rear axle:
% 
% Track width              m               1.715
% Toe Angle                deg             0.115          0.115         0.115
% Camber Angle             deg            -1.490         -1.490        -1.490
% 
% Normal force             N           10234.057       5117.028      5117.029
% Loaded tire radius       m               0.297          0.297         0.297
% Tire rate                N/mm          236.700        236.700       236.700
% ______________________________________________________________________________

% Initial values for road vehicle
         g = 9.81;
         Veh.ms  = 1841.860;            % sprung mass of the car in kg 
         Veh.mw = 266.136;              % unsprung mass of the car in kg
         Veh.Jw = 2.592;                % inertia of wheel in kgm^2
         Veh.L = 2.971;                 % Wheelbase
         Veh.lambda = 0.4958;           % Ratio of CoG
         Veh.Lf = Veh.L*Veh.lambda;     % Distance from front axis to CoG
         Veh.Lr = Veh.L*(1-Veh.lambda); % Distance from front axis to CoG
         Veh.h = 0.578;                 % CoG height from ground
         Veh.r = 0.297;                 % rolling radius for tire in m
         
        Veh.mu = [1 0.7 0.3];           % friction coefficients [dry, rain and snow]
        K_em = [230 160 80];
        K_brake = [600 500 400];
        tire_leg={'\mu = 1','\mu = 0.7','\mu = 0.3'};
        switch mu_select
            case {1}
                K_em = K_em(1)
                K_brake = K_brake(1)
                Veh.mu = Veh.mu(1)
            case {2}
                K_em = K_em(2)
                K_brake = K_brake(2)
                Veh.mu = Veh.mu(2)
            case {3}
                K_em = K_em(3)
                K_brake = K_brake(3)
                Veh.mu = Veh.mu(3)
        end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here you can ad your control parameters for the electric motor

Kp_em = 0;
Ki_em = 0;
Kd_em = 0;

% Here you can ad your control parameters for the brake 

Kp_brake = 0;
Ki_brake = 0;
Kd_brake = 0;

% Run your model
out = sim('SD2231_Lab1_2022.slx')

% Calculate the X-Y data, plot and estimate the optimum slip value
% Plot longitudinal slip on X-axis and used friction mu on y-axis

%% Task 2c
close all
mu_fl = -out.FxFL./out.FzFL;
mu_fr = -out.FxFR./out.FzFR;
mu_rl = -out.FxRL./out.FzRL;
mu_rr = -out.FxRR./out.FzRR;


FL_slip = out.FL_slip;
FR_slip = out.FR_slip;
RL_slip = out.RL_slip;
RR_slip = out.RR_slip;


figure
subplot(2,2,1)
plot(FL_slip(7660:8000),mu_fl(7660:8000))
grid on; title('FL')
subplot(2,2,2)
plot(FR_slip(7660:8000),mu_fr(7660:8000))
grid on; title('FR')
subplot(2,2,3)
plot(RL_slip(7660:8000),mu_rl(7660:8000))
grid on; title('RL')
subplot(2,2,4)
plot(RR_slip(7660:end),mu_rr(7660:end))
grid on; title('RR')


figure
plot(RR_slip(7660:8156),mu_rr(7660:8156))
grid on; title('RR')
%% Figures for report
figure
plot(RR_slip(175:3630),-mu_rr(175:3630),'LineWidth',1)
grid on
title('Traction')
tmp = mu_rr(175:3630);
[m,idx] = max(abs(tmp));
xlabel('Longitudinal slip $s$','Interpreter','latex')
ylabel('Friction coefficient $\mu$','Interpreter','latex')
hold on
plot(RR_slip(175+idx),-mu_rr(175+idx),'r*')
legend('RR','Max RR')
dim = [0.2 0.5 0.2 0.2];

str = {sprintf('Max point (%f,%f)', RR_slip(175+idx),-mu_rr(175+idx))};
annotation('textbox',dim,'String',str,'FitBoxToText','on');


figure
plot(RR_slip(7660:8156),mu_rr(7660:8156),'LineWidth',1)
grid on
title('Braking')
tmp = mu_rr(7660:8156);
[m,idx] = max(abs(tmp));
legend('RR')
xlabel('Longitudinal slip $s$','Interpreter','latex')
ylabel('Friction coefficient $\mu$','Interpreter','latex')
hold on
plot(RR_slip(7660+idx),mu_rr(7660+idx),'r*')
dim = [0.2 0.5 0.42 0.42];

str = {sprintf('Max point (%f,%f)', RR_slip(7660+idx),mu_rr(7660+idx))};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

hold on
plot(FL_slip(7660:8000),mu_fl(7660:8000),'LineWidth',1)
grid on
tmp = mu_fl(7660:8000);
[m,idx] = max(abs(tmp));

xlabel('Longitudinal slip $s$','Interpreter','latex')
ylabel('Friction coefficient $\mu$','Interpreter','latex')
hold on
plot(FL_slip(7660+idx),mu_fl(7660+idx),'b*')
legend('RR','Max RR','FL','Max FL')
dim = [0.2 0.5 0.17 0.17];

str = {sprintf('Max point (%f,%f)', FL_slip(7660+idx),mu_fl(7660+idx))};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

%% Ziegler-Nichols open-loop step response
%plant is fed with a step response on Throttle/Brake in open loop

%TRACTION
clf
figure
slip_rl = out.RL_slip(2050:7000);
time = out.tout(2050:7000);
plot(time,slip_rl,'LineWidth',1)
grid on
d2_rl = diff(diff(out.RL_slip(2050:7002)))/dt;
%plot(time,d2)
hold on

ds_rl=diff(slip_rl)./diff(time);
%I = find(time == 2.085);
d1_rl = diff(out.RL_slip(2050:7002))/dt;
k_rl=10; % point number 10
tang_rl=(time-time(k_rl))*d1_rl(k_rl)+slip_rl(k_rl);
plot(time,tang_rl,'LineWidth',1)
yline(max(out.RL_slip(2050:7002)),'--','LineWidth',1)
xline(2.0475,'--'); xline(2.063,'--')
legend("RL slip","tangent")

ylim([-0.5 1.5]); xlim([2.045 2.07])
xlabel("Time [s]"); ylabel("Slip")
L_new = 2.0475 - 2;
T_new = 2.063-2.0475;

str = {sprintf('L =%f', L_new)};
dim = [0.2 0.5 0.17 0.17];
annotation('textbox',dim,'String',str,'FitBoxToText','on');

str = {sprintf('T =%f', T_new)};
dim2 = [0.7 0.5 0.17 0.17];
annotation('textbox',dim2,'String',str,'FitBoxToText','on');

%%

clf
figure
slip_rl_b = out.RL_slip(7500:end-2);
time = out.tout(7500:end-2);
plot(time,slip_rl_b,'LineWidth',1)

grid on
d2_rl_b= diff(diff(out.RL_slip(7500:end)))/dt;
%plot(time,d2)
hold on

ds_rl_b=diff(slip_rl_b)./diff(time);
%I = find(time == 2.085);
d1_rl_b = diff(out.RL_slip(7500:end))/dt;
k_rl_b=720; % point number 10
tang_rl_b=(time-time(k_rl_b))*d1_rl_b(k_rl_b)+slip_rl_b(k_rl_b);
plot(time,tang_rl_b,'LineWidth',1)
yline(max(out.RL_slip(7500:end)),'--','LineWidth',1)
xline(8.048,'--'); xline(8.335,'--')
legend("RL slip","tangent")

ylim([-0.5 1.5]); xlim([7.95 8.5])
xlabel("Time [s]"); ylabel("Slip")
L_new_rl_b = 8.048 - 8;
T_new_rl_b = 8.335-8.048;

str = {sprintf('L =%f', L_new_rl_b)};
dim = [0.2 0.5 0.17 0.17];
annotation('textbox',dim,'String',str,'FitBoxToText','on');

str = {sprintf('T =%f', T_new_rl_b)};
dim2 = [0.7 0.5 0.17 0.17];
annotation('textbox',dim2,'String',str,'FitBoxToText','on');
%% ZN braking FL (not working)

figure
slip_fl_b = out.FL_slip(7500:end-2);
time = out.tout(7500:end-2);
plot(time,slip_fl_b,'LineWidth',1)


grid on
d2_fl_b= diff(diff(out.FL_slip(7500:end)))/dt;
%plot(time,d2_fl_b)
hold on

ds_fl_b=diff(slip_fl_b)./diff(time);
%I = find(time == 2.085);
d1_fl_b = diff(out.FL_slip(7500:end))/dt;
k_fl_b=564; % point number 10
tang_fl_b=(time-time(k_fl_b))*d1_fl_b(k_fl_b)+slip_fl_b(k_fl_b);
plot(time,tang_fl_b,'LineWidth',1)
yline(max(out.FL_slip(7500:end)),'--','LineWidth',1)
xline(7.987,'--'); xline(8.112,'--')
legend("RL slip","tangent")

ylim([-0.5 1.5]); xlim([7.95 8.15])
xlabel("Time [s]"); ylabel("Slip")
L_new_fl_b = 7.987 - 8;
T_new_fl_b = 8.112-7.987;

str = {sprintf('L =%f', L_new_fl_b)};
dim = [0.2 0.5 0.17 0.17];
annotation('textbox',dim,'String',str,'FitBoxToText','on');

str = {sprintf('T =%f', T_new_fl_b)};
dim2 = [0.7 0.5 0.17 0.17];
annotation('textbox',dim2,'String',str,'FitBoxToText','on');
