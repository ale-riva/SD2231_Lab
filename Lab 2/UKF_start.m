%----------------------------------------------------------------
% Template created for the course SD2231 by Mikael Nybacka 2019
% Following file is the start file for the state estimation using
% Uncented Kalman Filter (UKF).
%----------------------------------------------------------------
clear all;
close all;
clc;
addpath('scripts')
addpath('logged_data')
disp(' ');

% Set global variables so that they can be accessed from other matlab
% functions and files
global lf lr Ly_relax Cf Cr mass Mu Iz vbox_file_name g L

%----------------------------
% LOAD DATA FROM VBOX SYSTEM
%----------------------------
%Test for SD2231 20180912

%Cloudy 18 degrees celsius
%2.8 bar in all four tyres
%Two persons in the front

%vbox_file_name='S90__035.VBO';   %Standstill

%Need smoothing to get good plot
vbox_file_name='S90__036.VBO';   %Circular driving to the left, radius=8m
%Good estimation
%vbox_file_name='S90__038.VBO';  %Slalom, v=30km/h
% Need smooting but will not go down all the way
% vbox_file_name='S90__040.VBO';  %Step steer to the left, v=100km/h
%  Too large amplitude, need increasing Cf Cr, need negative ay
% vbox_file_name='S90__041.VBO';  %Frequency sweep, v=50km/h


vboload
%  Channel 1  = satellites
%  Channel 2  = time
%  Channel 3  = latitude
%  Channel 4  = longitude
%  Channel 5  = velocity kmh
%  Channel 6  = heading
%  Channel 7  = height
%  Channel 8  = vertical velocity kmh
%  Channel 9  = long accel g
%  Channel 10 = lat accel g
%  Channel 11 = glonass_sats
%  Channel 12 = gps_sats
%  Channel 13 = imu kalman filter status
%  Channel 14 = solution type
%  Channel 15 = velocity quality
%  Channel 16 = event 1 time
%  Channel 17 = latitude_raw
%  Channel 18 = longitude_raw
%  Channel 19 = speed_raw
%  Channel 20 = heading_raw
%  Channel 21 = height_raw
%  Channel 22 = vertical_velocity_raw
%  Channel 23 = rms_hpos
%  Channel 24 = rms_vpos
%  Channel 25 = rms_hvel
%  Channel 26 = rms_vvel
%  Channel 27 = poscov_xx
%  Channel 28 = poscov_yy
%  Channel 29 = poscov_zz
%  Channel 30 = velcov_xx
%  Channel 31 = velcov_yy
%  Channel 32 = velcov_zz
%  Channel 33 = t1
%  Channel 34 = accrpedlln 
%  Channel 35 = engn
%  Channel 36 = ptgearact
%  Channel 37 = pttqatw_fl 
%  Channel 38 = pttqatw_fr 
%  Channel 39 = swa
%  Channel 40 = brkpedlpsd 
%  Channel 41 = vehspdlgt
%  Channel 42 = flwhlspd
%  Channel 43 = frwhlspd
%  Channel 44 = rlwhlspd
%  Channel 45 = rrwhlspd
%  Channel 46 = algt1
%  Channel 47 = alat1
%  Channel 48 = rollrate1
%  Channel 49 = yawrate1
%  Channel 50 = true_head
%  Channel 51 = slip_angle 
%  Channel 52 = lat._vel.
%  Channel 53 = roll_angle 
%  Channel 54 = lng._vel.
%  Channel 55 = slip_cog
%  Channel 56 = yawrate
%  Channel 57 = x_accel
%  Channel 58 = y_accel
%  Channel 59 = temp
%  Channel 60 = pitchrate
%  Channel 61 = rollrate
%  Channel 62 = z_accel
%  Channel 63 = roll_imu
%  Channel 64 = pitch_ang. 
%  Channel 65 = yaw_rate
%  Channel 66 = slip_fl
%  Channel 67 = slip_fr
%  Channel 68 = slip_rl
%  Channel 69 = slip_rr
%  Channel 70 = true_head2 
%  Channel 71 = head_imu
%  Channel 72 = pitch_imu
%  Channel 73 = pos.qual.
%  Channel 74 = lng_jerk
%  Channel 75 = lat_jerk
%  Channel 76 = head_imu2

%%
%-----------------------------------
% SET VEHICLE DATA FOR THE VOLVO S90
%-----------------------------------
Rt=0.35;            % Tyre radius (m)
L=2.941;            % Wheel base (m)
lf=1.65;            % Distance from CoG to front axis (m)
lr=L-lf;            % Distance from CoG to rear axis (m)
mass=2010.5;        % Mass (kg)
Iz=3089;            % Yaw inertia (kg-m2)
tw=1.617;           % Track width (m)
h_cog = 0.570;      % Height of CoG above ground
Ratio=16.3;         % Steering gear ratio
Cf=200000;          % Lateral stiffness front axle (N)
Cr=250000;          % Lateral stiffness rear axle (N)
Lx_relax=0.05;      % Longitudinal relaxation lenth of tyre (m)
Ly_relax=0.15;      % Lateral relaxation lenth of tyre (m)
Roll_res=0.01;      % Rolling resistance of tyre
rollGrad=4.5;       % rollangle deg per latacc
rx=-0.29;           % distance from CoG to IMU x-axle
ry=0;               % distance from CoG to IMU y-axle
rz=0;               % distance from CoG to IMU z-axle

Ts = 0.01;

%--------------------------------------
% SET ENVIRONEMNTAL PARAMETERS FOR TEST
%--------------------------------------
Mu=0.85;            % Coefficient of friction
g=9.81;             % Gravity constant

%--------------------------------------------
% SET VARIABLES DATA FROM DATA READ FROM FILE
%--------------------------------------------

Time            = vbo.channels(1, 2).data-vbo.channels(1, 2).data(1,1);
yawRate_VBOX    = vbo.channels(1, 56).data.*(-pi/180);%VBOX z-axis is pointing downwards, hence (-)
vx_VBOX         = vbo.channels(1, 54).data./3.6;
vy_VBOX         = vbo.channels(1, 52).data./3.6;
SteerAngle      = vbo.channels(1, 39).data./Ratio;
ax_VBOX         = vbo.channels(1, 57).data.*g;
ay_VBOX         = vbo.channels(1, 58).data.*g;
Beta_VBOX       = (vy_VBOX + rx*yawRate_VBOX)./vx_VBOX;

% Time_cut = Time;
Time_cut = Time(1001:end-400);
yawRate_VBOX    = yawRate_VBOX(1001:end-400);
vx_VBOX = vx_VBOX(1001:end-400);
vy_VBOX = vy_VBOX(1001:end-400);
SteerAngle = SteerAngle(1001:end-400);
ax_VBOX = ax_VBOX(1001:end-400);
ay_VBOX = ay_VBOX(1001:end-400);
Beta_VBOX = Beta_VBOX(1001:end-400);
slip_rr        = vbo.channels(1, 69).data;
rrwhlsp        = vbo.channels(1,45).data;
slip_rr_formula =-(Rt.*rrwhlsp(size(vx_VBOX,2),1)' - vx_VBOX)./Rt.*rrwhlsp(size(vx_VBOX,2),1)';
figure
plot(slip_rr_formula)

n = length(Time_cut);
dt = Time(2)-Time(1);
%%
%----------------------------------------------
% SET MEASUREMENT AND PROCESS NOICE COVARIANCES
%----------------------------------------------
% Use as starting value 0.1 for each of the states in Q matrix
% Q=diag([0.1;0.1;0.1]);    % Default tuning
%Q=diag([0.04;0.07;0.03]); % Better result, from std of sim1. Final tune

% For extended system
% Q=diag([2;6;0.27]);

model = 1;
%0 for linear tyre model
%1 for brush tyre model
if model == 0
%     Q=diag([0.04;0.07;0.03]); 
%     R=diag([2.8;0.2;0.08]); 
      Q=diag([0.1;0.1;0.1]);
      R=diag([0.01,0.01,0.01])
elseif model == 1
%     Q=diag([2;6;0.27]);
%     R=diag([6.88;3.53;0.15]);
      Q=diag([0.04;0.01;0.03]); 
      R=diag([2.8;0.2;0.08]);
      Cf = 160000;
      Cr = 200000;
end

% Use as starting value 0.01 for each of the measurements in R matrix
% R=diag([0.01,0.01,0.01]);   % Default tuning
%R=diag([2.8;0.2;0.08]); % Better result, from std of sim1. Final tune

% For extended system
% R=diag([6.88;3.53;0.15]);


%--------------------------------------------------
% SET INITIAL STATE AND STATE ESTIMATION COVARIANCE
%--------------------------------------------------
x_0= [0.01;0;0];  %Don't set vx to zero, otherwise divison-by-zero will occur
P_0= diag([0.01;0.1;0.01]);

% For extended system
% vx_0 = 0.01;    % Don't set to zero
% vy_0 = 0;
% yawrate_0 = 0;
% x_0= [vx_0;vy_0;yawrate_0];
% P_0= diag([0.1;0.1;0.1]);


%-----------------------
% INITIALISING VARIABLES
%-----------------------


%Parameters that might be needed in the measurement and state functions are added to predictParam
predictParam.dt=dt; 

% Handles to state and measurement model functions.
if model == 0
    state_func_UKF = @Vehicle_state_eq;
    meas_func_UKF = @Vehicle_measure_eq;
elseif model == 1
    state_func_UKF = @Vehicle_state_eq_new2;
    meas_func_UKF = @Vehicle_measure_eq_new2;
end

%-----------------------
% FILTERING LOOP FOR UKF 
%-----------------------
disp(' ');
disp('Filtering the signal with UKF...');

x=x_0;
P = P_0;
x_mat = x;

% smooth(ay_VBOX(i)+rx*(yawRate_VBOX(i)-yawRate_VBOX(i-1))/dt);
for i = 2:n
    [x,P] = ukf_predict1(x,P,state_func_UKF,Q,SteerAngle(i),0.8,2,0); %Alpha,Beta,Kappa 0.8, 2, 0

    meas_Y = [smooth(vx_VBOX(i));smooth(ay_VBOX(i)+rx*(yawRate_VBOX(i)-yawRate_VBOX(i-1))/dt);smooth(yawRate_VBOX(i))];
    [x,P] = ukf_update1(x,P,meas_Y,meas_func_UKF,R,SteerAngle(i),0.8,2,0);
    x_mat = [x_mat,x];

    % ad your predict and update functions, see the scripts ukf_predict1.m
    % and ukf_update1.m

    
    if i==round(n/4)
        disp(' ');
        disp('1/4 of the filtering done...');
        disp(' ');
    end
    if i==round(n/2)
        disp(' ');
        disp('1/2 of the filtering done...');
        disp(' ');
    end
    if i==round(n*(3/4))
        disp(' ');
        disp('3/4 of the filtering done... Stay tuned for the results...');
        disp(' ');
    end
end
disp('Filtering done');
%%
%----------------------------------------
% CALCULATE THE SLIP ANGLE OF THE VEHICLE
%----------------------------------------

Beta_vy = atan2(x_mat(2,:),x_mat(1,:));
%---------------------------------------------------------
% CALCULATE THE ERROR VALES FOR THE ESTIMATE OF SLIP ANGLE
%---------------------------------------------------------
% Beta_VBOX_smooth=smooth(Beta_VBOX(1400:4728),0.01,'rlowess'); 
% [e_beta_mean,e_beta_max,time_at_max,error] = errorCalc(Beta_vy(1700:2152),Beta_VBOX_smooth);
% disp(' ');
% fprintf('The MSE of Beta estimation is: %d \n',e_beta_mean);
% fprintf('The Max error of Beta estimation is: %d \n',max(e_beta_max));

%-----------------
% PLOT THE RESULTS
%-----------------
% clf
% figure(1)
% subplot(2,2,1)
% plot(Time,vx_VBOX)
% hold on
% plot(Time,x_mat(1,:))
% title("vx velocity")
% subplot(2,2,2)
% plot(Time,smooth(vy_VBOX))
% hold on
% plot(Time,smooth(x_mat(2,:)))
% title("vy velocity")
% ylim([-1,1])
% subplot(2,2,3)
% plot(Time,yawRate_VBOX)
% hold on
% plot(Time,x_mat(3,:))
% title("yawrate")
% subplot(2,2,4)
% plot(Time,smooth(Beta_VBOX))
% hold on
% plot(Time,smooth(Beta_vy))
% title("Side slip")
% ylim([-0.1,0.3])
% disp("Plotting coming")

% figure(2)
% plot(Time,Beta_VBOX,'LineWidth',1)
% hold on
% plot(Time,smooth(Beta_vy),'LineWidth',1)
% hold on
% plot(Time(1200:end-200),out.Betay_mod.Data,'LineWidth',1)
% hold on
% plot(Time(1200:end-200),out.Betay_kin.Data,'LineWidth',1)
% hold on
% plot(Time(1200:end-200),out.Betay_wf.Data,'LineWidth',1)
% hold on
% plot(Time(1200:end-200),out.Betay_wf_var.Data,'LineWidth',1,'Color','g')
% title("Side slip, sim 2")
% legend("VBOX","UKF", "Model","Kinetic","WF, fixed T","WF, var T")
% % legend("VBOX","UKF","WF, var T")
% xlim([Time(1200),Time(end-200)])
% ylim([-0.1,0.5])
% xlabel("Time")
% ylabel("Slip amplitude")
% grid on

%%
vy_COG = vy_VBOX  + rx.*yawRate_VBOX;

figure
plot(Time_cut,vx_VBOX)
hold on
plot(Time_cut,x_mat(1,:))
legend("meas","UKF")
title("vx velocity")
figure
plot(Time_cut,smooth(vy_COG))
hold on
plot(Time_cut,smooth(x_mat(2,:)))
title("vy velocity")
legend("meas","UKF")
%ylim([-1,1])
figure
plot(Time_cut,yawRate_VBOX)
hold on
plot(Time_cut,x_mat(3,:))
title("yawrate")
legend("meas","UKF")
figure
plot(Time_cut,smooth(Beta_VBOX))
hold on
plot(Time_cut,smooth(Beta_vy))
title("Side slip")
legend("meas","UKF")
%ylim([-0.1,0.3])
disp("Plotting coming")

if model == 0
    vx_lin = x_mat(1,:);
    vy_lin = x_mat(2,:);
    yawRate_lin = x_mat(3,:);
    beta_lin = Beta_vy;
elseif model ==1
    vx_ext = x_mat(1,:);
    vy_ext = x_mat(2,:);
    yawRate_ext = x_mat(3,:);
    beta_ext = Beta_vy;
end

%% MSE and MAX
start_time_index = 1001;
time_start = Time(start_time_index);


sim1 = [4 37.28];
sim2 = [5.5 21.78];
sim3 = [7 11.510];
sim4 = [0 27.23];

switch vbox_file_name
    case 'S90__036.VBO'
        sim_txt = "sim1";
        sim_ok = sim1;
    case 'S90__038.VBO'
        sim_txt = "sim2";
        sim_ok = sim2;
    case 'S90__040.VBO'
        sim_txt = "sim3";
        sim_ok = sim3;
    case 'S90__041.VBO'
        sim_txt = "sim4";
        sim_ok = sim4;
end
sim_ok
sim = sim_ok.*100+1;
sim_truth = (sim_ok +time_start).*100+1;

if model == 0
    vx_ukf = vx_lin;
    vy_ukf = vy_lin; 
    yawRate_ukf = yawRate_lin; 
    beta_ukf = beta_lin ;
elseif model ==1
    vx_ukf = vx_ext;
    vy_ukf = vy_ext; 
    yawRate_ukf = yawRate_ext; 
    beta_ukf = beta_ext ;
end

figure
plot(smooth(beta_ukf(sim(1):sim(2))));
hold on
plot(smooth(Beta_VBOX(sim(1):sim(2))));
title("check")
%Beta_VBOX_smooth(start_time_index:(start_time_index+length(out.Betay_mod.Data)-1),2)
[error_mean_vx,error_max_vx,time_at_max_vx,~] = errorCalc(vx_ukf(sim(1):sim(2)),vx_VBOX(sim(1):sim(2))');
fprintf("VX: MSE=%e; MAX=%e\n",error_mean_vx,error_max_vx)

[error_mean_vy,error_max_vy,time_at_max_vy,~] = errorCalc(smooth(vy_ukf(sim(1):sim(2))),smooth(vy_COG(sim(1):sim(2))'));
fprintf("VY: MSE=%e; MAX=%e\n",error_mean_vy,error_max_vy)

[error_mean_yr,error_max_yr,time_at_max_yr,~] = errorCalc(yawRate_ukf(sim(1):sim(2)),yawRate_VBOX(sim(1):sim(2))');
fprintf("YR: MSE=%e; MAX=%e\n",error_mean_yr,error_max_yr)

[error_mean_beta,error_max_beta,time_at_max_beta,~] = errorCalc(smooth(beta_ukf(sim(1):sim(2))),smooth(Beta_VBOX(sim(1):sim(2))'));
fprintf("BETA: MSE=%e; MAX=%e\n",error_mean_beta,error_max_beta)

%%

fileID = fopen('task2.txt','a');
fprintf(fileID,"\n\n simulation: %s\n",sim_txt)
if model == 0
    fprintf(fileID,"MODEL linear\n");
elseif model ==1
    fprintf(fileID,"MODEL extended\n");
end
fprintf(fileID,"R:[%f,%f,%f]\n",diag(R).');
fprintf(fileID,"Cf:%d,   Cr:%d\n",Cf,Cr);
%fprintf(fileID,"time interval: [%f %f]\n",start_time_index)
fprintf(fileID,"VX: MSE=%e; MAX=%e\n",error_mean_vx,error_max_vx)
fprintf(fileID,"VY: MSE=%e; MAX=%e\n",error_mean_vy,error_max_vy)
fprintf(fileID,"YR: MSE=%e; MAX=%e\n",error_mean_yr,error_max_yr)
fprintf(fileID,"BETA: MSE=%e; MAX=%e\n",error_mean_beta,error_max_beta)
fclose(fileID);
    