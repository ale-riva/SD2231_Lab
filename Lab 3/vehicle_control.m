clear all
close all
clc


%% Task 8.1 
w=logspace(-2,3,2*100040);

m = 22000; J = 700000;
c = 40000;
k = 600000;
k1 = k; k2= k;
c1 = c; c2 = c;
L = 6;  L1= L; L2= L;
M_vh = [m 0;0 J];
K_vh = [(k1+k2) (k2*L2-k1*L1);
        (k2*L2-k1*L1) (k2*L2^2+k1*L1^2)];
C_vh = [(c1+c2) (c2*L2-c1*L1); (c2*L2-c1*L1) (c2*L2^2+c1*L1^2)];

wn_und_vh = sqrt(eig(K_vh,M_vh))

A_ss_2dof_vh = [0  1 0 0;
    -(k1+k2)/m ,-(c1+c2)/m,(k1*L1-k2*L2)/m,(c1*L1-c2*L2)/m;
     0, 0, 0 ,1;
     (k1*L1-k2*L2)/J, (c1*L1-c2*L2)/J,(-k1*L1^2-k2*L2^2)/J,(-c1*L1^2-c2*L2^2)/J];
eig(A_ss_2dof_vh)
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

%% check with bode polot for z, chi from input z1

veh_2dof_tf = tf(ss(A_ss_2dof_vh,B_ss_2dof_vh,C_ss_2dof_vh,D_ss_2dof_vh));
veh_2dof_tf.InputName = ["z_{w1}"," z_{w1} dot","z_{w2}"," z_{w2} dot"];
veh_2dof_tf.OutputName = ["zeta","chi"];
[A_veh_2dof,phi_veh_2dof]=bode(veh_2dof_tf,w); 

[wn_veh,r_veh] = damp(veh_2dof_tf)
tmp1 = reshape(A_veh_2dof(1,1,:),[],size(A_veh_2dof(1,1,:), 3));
tmp2 = reshape(A_veh_2dof(2,1,:),[],size(A_veh_2dof(2,1,:), 3));

figure
loglog(w,tmp1)
hold on
loglog(w,tmp2)
grid on
hold on
xline(1*2*pi,"r--")
hold on
xline(8*2*pi,"p--")
legend("z zw1","chi zw1","1 Hz","8 hz");



%% Task 9.1
A_veh_sh = [0 1 0 0;
            -(k2+k1)/m, 0, (k1*L1-k2*L2)/m, 0;
            0  0  0  1;
            (k1*L1-k2*L2)/J, 0, -(k2*L2^2+k1*L1^2)/J, 0];
B_veh_sh = [ 0 0 0 0;
             k1/m k2/m -1/m -1/m;
             0      0    0   0;
             -k1*L1/J k2*L2/J L1/J    -L2/J];
C_veh_sh = [0 1 0 0;0 0 0 1];
D_veh_sh = [0 0 0 0;0 0 0 0];

cz=15e4;    %Lower force for impulse
c_chi=8e5;

% Close loop check
A_veh_sh_cl = [0   1    0     0;
               -(k1+k2)/m  -cz/m   (k1*L1-k2*L2)/m  0;
               0           0           0           1;
               (k1*L1-k2*L2)/J   0     -(k1*L1^2+k2*L2^2)/J  -c_chi/J];
B_veh_sh_cl = [0 0; k1/m k2/m;0 0 ;-k1*L1/J k2*L2/J];
C_veh_sh_cl = [1 0 0 0; 0 0 1 0];
D_veh_sh_cl = [0 0; 0 0];

veh_sh_tf_cl = tf(ss(A_veh_sh_cl,B_veh_sh_cl,C_veh_sh_cl,D_veh_sh_cl));
[mag_veh_sh_cl,phi_veh_sh_cl] =bode(veh_sh_tf_cl,w);


%bode plot for input zw1 of the CLOSED LOOP 
tmp3 = reshape(mag_veh_sh_cl(1,1,:),[],size(A_veh_2dof(2,1,:), 3));
tmp4 = reshape(mag_veh_sh_cl(2,1,:),[],size(A_veh_2dof(2,1,:), 3));

% figure
% subplot(2,1,1)
% loglog(w,tmp1)
% hold on
% loglog(w,tmp3)
% hold on
% xline(1*2*pi,"r--")
% hold on
% xline(8*2*pi,"p--")
% title("z")
% legend("ol","cl","1 hz","8 hz")
% grid on
% subplot(2,1,2)
% loglog(w,tmp2)
% hold on
% loglog(w,tmp4)
% hold on
% xline(1*2*pi,"r--")
% hold on
% xline(8*2*pi,"p--")
% grid on
% title("chi")
% legend("ol","cl","1 hz","8 hz")


% 
% run("SS_task9.m");
% out9_1 = out9

%% Task 9.3
out9_1 = sim("ale_SS_task9.slx",'StartTime','0','StopTime','20','FixedStep','0.01');
time = out8_2.tout;
exc1_z_resp_veh_sh = out9_1.z_exc1;
exc1_chi_resp_veh_sh = out9_1.chi_exc1;
exc2_8_z_resp_veh_sh = out9_1.z_exc2_8;
exc2_8_chi_resp_veh_sh = out9_1.chi_exc2_8;
exc2_1_z_resp_veh_sh = out9_1.z_exc2_1;
exc2_1_chi_resp_veh_sh = out9_1.chi_exc2_1;

figure
subplot(2,1,1)
plot(time,exc1_z_resp_veh_sh,'r', 'LineWidth',1)
hold on
plot(time,exc1_z_resp_2dof_vh,'r--')
title("z")
legend("SH","Passive");
xlim([0,10]);
grid on
subplot(2,1,2)
plot(time,exc1_chi_resp_veh_sh,'b', 'LineWidth',1)
hold on
plot(time,exc1_chi_resp_2dof_vh,'b--')
sgtitle("response of vehicle to EXC1 - SH")
title("$\chi$",'Interpreter','latex')
legend("SH","Passive");         
xlim([0,10]);
grid on

figure
subplot(2,1,1)
plot(time,exc2_1_z_resp_veh_sh,'r', 'LineWidth',1)
hold on
plot(time,exc2_1_z_resp_2dof_vh,'r--')
title("z")
legend("SH","Passive");
xlim([0,10]);
grid on
subplot(2,1,2)
plot(time,exc2_1_chi_resp_veh_sh,'b', 'LineWidth',1)
hold on
plot(time,exc2_1_chi_resp_2dof_vh,'b--')
sgtitle("response of vehicle to EXC2 (1Hz) - SH")
title("$\chi$",'Interpreter','latex')
legend("SH","Passive");         
xlim([0,10]);
grid on

figure
subplot(2,1,1)
plot(time,exc2_8_z_resp_veh_sh,'r', 'LineWidth',1)
hold on
plot(time,exc2_8_z_resp_2dof_vh,'r--')
title("z")
legend("SH","Passive");
xlim([0,10]);
grid on
subplot(2,1,2)
plot(time,exc2_8_chi_resp_veh_sh,'b', 'LineWidth',1)
hold on
plot(time,exc2_8_chi_resp_2dof_vh,'b--')
sgtitle("response of vehicle to EXC2 (8Hz) - SH")
title("$\chi$",'Interpreter','latex')
legend("SH","Passive");         
xlim([0,10]);
grid on

