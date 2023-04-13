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
Kp_br_f = 10*T/L;
Ki_br_f = 10*T/(L^2);
Kd_br_f = 0.1 *T;

%final values morning 13/04
Kp_br_r = 2*T/L;
Ki_br_r = 2*T/(L^2);
Kd_br_r = 2*T;


%final values afternoon 13/04
% Kp_br_r = 2*T/L
% Ki_br_r = 0.7*T/(L^2)
% Kd_br_r = 4*T;



Kp_tr =1*T/L;
Ki_tr =0.4*T/(L^2);
Kd_tr = 2*T;

N = 100;
delay = 250;


K_awup = Ki_br_r/Kp_br_r;
ref_tr = 0.06;
ref_br_f = 0.15;
ref_br_r = 0.17; 
threshold = 0.01;
threshold_tc = ref_tr;
threshold_br_f = ref_br_f;
threshold_br_r = ref_br_r-0.05;
step_time = 2;
vx_threshold = 0;
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

%%
clf
figure
slip = out.RL_slip(2050:7000);
time = out.tout(2050:7000);
plot(time,slip)
hold on
grid on
d2 = diff(diff(out.RL_slip(2050:7002)))/dt;
plot(time,d2)
hold on


ds=diff(slip)./diff(time);
I = find(time == 2.085);
d1 = diff(out.RL_slip(2050:7002))/dt;
k=I; % point number 220
tang=(time-time(k))*d1(k)+slip(k);
plot(time,tang)
yline(max(out.RL_slip(2050:7002)))
legend("slip","2d","t")
%scatter(time(k),slip(k))


