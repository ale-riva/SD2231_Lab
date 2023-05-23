clear all
close all
clc


%% Task 8.1 
w=logspace(-2,3,2*100040);

m = 22000; J = 700000;
c = 40000;
k = 600000;
k1 = k; k2= k;
c1 = c,c2 = c;
L = 6;  L1= L; L2= L;
M_vh = [m 0;0 J];
K_vh = [(k1+k2) (k2*L2-k1*L1);
        (k2*L2-k1*L1) (k2*L2^2+k1*L1^2)];
wn_und_vh = sqrt(eig(K_vh,M_vh))

A_ss_2dof_vh = [0  1 0 0;
    -(k1+k2)/m ,-(c1+c2)/m,(k1*L1-k2*L2)/m,(c1*L1-c2*L2)/m;
     0, 0, 0 ,1;
     (k1*L1-k2*L2)/J, (c1*L1-c2*L2)/J,(-k1*L1^2-k2*L2^2)/J,(-c1*L1^2-c2*L2^2)/J];
B_ss_2dof_vh = [0 0 0 0;
                k1/m c1/m  k2/m  c2/m;
                 0     0     0     0;
                 -(k1*L1)/m  -c1*L1/J k2*L2/J c2*L2/J];
C_ss_2dof_vh = [1 0 0 0;
                0 0 1 0];
D_ss_2dof_vh = [0 0 0 0; 0 0 0 0];


%% Task 8.2
veh_2dof_tf = tf(ss(A_ss_2dof_vh,B_ss_2dof_vh,C_ss_2dof_vh,D_ss_2dof_vh));
veh_2dof_tf.InputName = ["z_{w1}"," z_{w1} dot","z_{w2}"," z_{w2} dot"];
veh_2dof_tf.OutputName = ["zeta","chi"];
out8_2 = sim("SS_task8_2.slx",'StartTime','0','StopTime','20','FixedStep','0.01');

time = out8_2.tout;
exc1_z_resp_2dof_vh = out8_2.z_exc1;
exc1_chi_resp_2dof_vh = out8_2.chi_exc1;
exc2_8_z_resp_2dof_vh = out8_2.z_exc2_8;
exc2_8_chi_resp_2dof_vh = out8_2.chi_exc2_8;
exc2_1_z_resp_2dof_vh = out8_2.z_exc2_1;
exc2_1_chi_resp_2dof_vh = out8_2.chi_exc2_1;

figure
plot(time,exc1_z_resp_2dof_vh)
hold on
plot(time,exc1_chi_resp_2dof_vh)
title("response of vehicle to EXC1 (impulse on zw1)")
legend("z","$\chi$",'Interpreter','latex')         
xlim([0,10]);
grid on

figure
plot(time,exc2_1_z_resp_2dof_vh)
hold on
plot(time,exc2_1_chi_resp_2dof_vh)
title("response of vehicle to EXC2 (1Hz)")
legend("z","$\chi$",'Interpreter','latex')         
xlim([0,10]);
grid on

figure
plot(time,exc2_8_z_resp_2dof_vh)
hold on
plot(time,exc2_8_chi_resp_2dof_vh)
title("response of vehicle to EXC2 (8Hz)")
legend("z","$\chi$",'Interpreter','latex')         
xlim([0,10]);
grid on

%%
veh_2dof_tf = tf(ss(A_ss_2dof_vh,B_ss_2dof_vh,C_ss_2dof_vh,D_ss_2dof_vh));
veh_2dof_tf.InputName = ["z_{w1}"," z_{w1} dot","z_{w2}"," z_{w2} dot"];
veh_2dof_tf.OutputName = ["zeta","chi"];
[A_veh_2dof,phi_veh_2dof]=bode(veh_2dof_tf,w); 
%%
[wn_veh,r_veh] = damp(veh_2dof_tf)
%mag_lin_z_zw1 = squeeze(A_veh_2dof_und(1,1,:)); % Convert magnitude to linear scale
%mag_lin_chi_zw1 = squeeze(A_veh_2dof_und(2,1,:));
%[~, peak_indices_z_zw1] = findpeaks(mag_lin_z_zw1); % Find peak indices
%[~, peak_indices_chi_zw1] = findpeaks(mag_lin_chi_zw1); % Find peak indices
%peak_indices=[peak_indices_z_zw1,peak_indices_chi_zw1]
%nat_freq = w(peak_indices) % Natural frequencies

tmp1 = reshape(A_veh_2dof(1,1,:),[],size(A_veh_2dof(1,1,:), 3));
tmp2 = reshape(A_veh_2dof(2,1,:),[],size(A_veh_2dof(2,1,:), 3));

figure
loglog(w,tmp1)
hold on
loglog(w,tmp2)
grid on

hold on
xline(1*2*pi)
hold on
xline(8*2*pi)
legend("z zw1","chi zw1","wn z","wn chi");

tmp3 = reshape(A_veh_2dof(1,3,:),[],size(A_veh_2dof(1,1,:), 3));
tmp4 = reshape(A_veh_2dof(2,3,:),[],size(A_veh_2dof(2,1,:), 3));

figure
loglog(w,tmp3)
hold on
loglog(w,tmp4)
grid on

hold on
xline(1*2*pi)
hold on
xline(8*2*pi)
legend("z zw2","chi zw2","wn z","wn chi");

