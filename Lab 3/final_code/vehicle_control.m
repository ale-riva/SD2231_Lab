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
z_pas_exc1 = out8_2.z_exc1;
chi_pas_exc1 = out8_2.chi_exc1;
z_pas_exc2_8 = out8_2.z_exc2_8;
chi_pas_exc2_8 = out8_2.chi_exc2_8;
z_pas_exc2_1 = out8_2.z_exc2_1;
chi_pas_exc2_1 = out8_2.chi_exc2_1;

figure
plot(time,z_pas_exc1)
hold on
plot(time,chi_pas_exc1)
title("response of vehicle to EXC1 (impulse on zw1)")
legend("z","$\chi$",'Interpreter','latex')         
xlim([0,10]);
grid on

figure
plot(time,z_pas_exc2_1)
hold on
plot(time,chi_pas_exc2_1)
title("response of vehicle to EXC2 (1Hz)")
legend("z","$\chi$",'Interpreter','latex')         
xlim([0,10]);
grid on

figure
plot(time,z_pas_exc2_8)
hold on
plot(time,chi_pas_exc2_8)
title("response of vehicle to EXC2 (8Hz)")
legend("z","$\chi$",'Interpreter','latex')         
xlim([0,10]);
grid on

%% check with bode polot for z, chi from input z1

% veh_2dof_tf = tf(ss(A_ss_2dof_vh,B_ss_2dof_vh,C_ss_2dof_vh,D_ss_2dof_vh));
% veh_2dof_tf.InputName = ["z_{w1}"," z_{w1} dot","z_{w2}"," z_{w2} dot"];
% veh_2dof_tf.OutputName = ["zeta","chi"];
% [A_veh_2dof,phi_veh_2dof]=bode(veh_2dof_tf,w); 
% 
% [wn_veh,r_veh] = damp(veh_2dof_tf)
% tmp1 = reshape(A_veh_2dof(1,1,:),[],size(A_veh_2dof(1,1,:), 3));
% tmp2 = reshape(A_veh_2dof(2,1,:),[],size(A_veh_2dof(2,1,:), 3));
% 
% figure
% loglog(w,tmp1)
% hold on
% loglog(w,tmp2)
% grid on
% hold on
% xline(1*2*pi,"r--")
% hold on
% xline(8*2*pi,"p--")
% legend("z zw1","chi zw1","1 Hz","8 hz");



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
% tmp3 = reshape(mag_veh_sh_cl(1,1,:),[],size(A_veh_2dof(2,1,:), 3));
% tmp4 = reshape(mag_veh_sh_cl(2,1,:),[],size(A_veh_2dof(2,1,:), 3));



%% Task 9.3
out9_1 = sim("SS_task9.slx",'StartTime','0','StopTime','20','FixedStep','0.01');
time = out8_2.tout;
z_sh_exc1 = out9_1.z_exc1;
chi_sh_exc1 = out9_1.chi_exc1;
z_sh_exc2_8 = out9_1.z_exc2_8;
chi_sh_exc2_8 = out9_1.chi_exc2_8;
z_sh_exc2_1 = out9_1.z_exc2_1;
chi_sh_exc2_1 = out9_1.chi_exc2_1;

figure
subplot(2,1,1)
plot(time,z_sh_exc1,'r', 'LineWidth',1)
hold on
plot(time,z_pas_exc1,'r--')
title("z")
legend("SH","Passive");
xlim([0,10]);
grid on
subplot(2,1,2)
plot(time,chi_sh_exc1,'b', 'LineWidth',1)
hold on
plot(time,chi_pas_exc1,'b--')
sgtitle("response of vehicle to EXC1 - SH")
title("$\chi$",'Interpreter','latex')
legend("SH","Passive");         
xlim([0,10]);
grid on

figure
subplot(2,1,1)
plot(time,z_sh_exc2_1,'r', 'LineWidth',1)
hold on
plot(time,z_pas_exc2_1,'r--')
title("z")
legend("SH","Passive");
xlim([0,10]);
grid on
subplot(2,1,2)
plot(time,chi_sh_exc2_1,'b', 'LineWidth',1)
hold on
plot(time,chi_pas_exc2_1,'b--')
sgtitle("response of vehicle to EXC2 (1Hz) - SH")
title("$\chi$",'Interpreter','latex')
legend("SH","Passive");         
xlim([0,10]);
grid on

figure
subplot(2,1,1)
plot(time,z_sh_exc2_8,'r', 'LineWidth',1)
hold on
plot(time,z_pas_exc2_8,'r--')
title("z")
legend("SH","Passive");
xlim([0,10]);
grid on
subplot(2,1,2)
plot(time,chi_sh_exc2_8,'b', 'LineWidth',1)
hold on
plot(time,chi_pas_exc2_8,'b--')
sgtitle("response of vehicle to EXC2 (8Hz) - SH")
title("$\chi$",'Interpreter','latex')
legend("SH","Passive");         
xlim([0,10]);
grid on

% check on Maximum Fa1,Fa2
Fa1_1 = out9_1.exc1_Fa1;
Fa2_1 = out9_1.exc1_Fa2;
Fa1_2_1 = out9_1.exc2_1_Fa1;
Fa2_2_1 = out9_1.exc2_1_Fa2;
Fa1_2_8 = out9_1.exc2_8_Fa1;
Fa2_2_8 = out9_1.exc2_8_Fa2;

F_check = [Fa1_1;Fa2_1;Fa1_2_1;Fa2_2_1;Fa1_2_8;Fa2_2_8];
check_magn_SH = max(F_check)

%% Task 10 - Hinf
run("H_inf_ctrl.m")
z_hinf_exc1 = out10.z_exc1;
chi_hinf_exc1 = out10.chi_exc1;
z_hinf_exc2_1 = out10.z_exc2_1;
chi_hinf_exc2_1 = out10.chi_exc2_1;
z_hinf_exc2_8 = out10.z_exc2_8;
chi_hinf_exc2_8 = out10.chi_exc2_8;
time_hinf = out10.tout;
%%
figure
subplot(2,1,1)
plot(time_hinf,z_hinf_exc1,'LineWidth',1)
hold on
plot(time,z_pas_exc1,'--')
grid on
title("$z$ Excitation 1",'Interpreter','latex')
legend("H inf","Passive")

subplot(2,1,2)
plot(time_hinf,chi_hinf_exc1,'LineWidth',1)
hold on
plot(time,chi_pas_exc1,'--')
grid on
title("$\chi$ Excitation 1",'Interpreter','latex')
legend("H inf","Passive")


figure
subplot(2,1,1)
plot(time_hinf,z_hinf_exc2_1,'LineWidth',1)
hold on
plot(time,z_pas_exc2_1,'--')
grid on
title("$z$ Excitation 2.1",'Interpreter','latex')
legend("H inf","Passive")

subplot(2,1,2)
plot(time_hinf,chi_hinf_exc2_1,'LineWidth',1)
hold on
plot(time,chi_pas_exc2_1,'--')
grid on
title("$\chi$ Excitation 2.1",'Interpreter','latex')
legend("H inf","Passive")

figure
subplot(2,1,1)
plot(time_hinf,z_hinf_exc2_8,'LineWidth',1)
hold on
plot(time,z_pas_exc2_8,'--')
grid on
title("$z$ Excitation 2.8",'Interpreter','latex')
legend("H inf","Passive")

subplot(2,1,2)
plot(time_hinf,chi_hinf_exc2_8,'LineWidth',1)
hold on
plot(time,chi_pas_exc2_8,'--')
grid on
title("$\chi$ Excitation 2.8",'Interpreter','latex')
legend("H inf","Passive")

Fa1_1_inf = out10.Fa1_exc1;
Fa2_1_inf = out10.Fa2_exc1;
Fa1_2_1_inf = out10.Fa1_exc2_1;
Fa2_2_1_inf = out10.Fa2_exc2_1;
Fa1_2_8_inf = out10.Fa1_exc2_8;
Fa2_2_8_inf = out10.Fa2_exc2_8;

F_check_inf = [Fa1_1_inf;Fa2_1_inf;Fa1_2_1_inf;Fa2_2_1_inf;Fa1_2_8_inf;Fa2_2_8_inf];
check_magn_inf = max(F_check_inf)


%% Task 11 - Summary
figure
plot(time_hinf,z_hinf_exc1,'LineWidth',1)
hold on
plot(time,z_sh_exc1,'LineWidth',1)
hold on
plot(time,z_pas_exc1,'m--')
grid on
title("$z$ Excitation 1",'Interpreter','latex')
legend("H inf","SH","Passive")
xlim([0,5])

figure
plot(time_hinf,chi_hinf_exc1,'LineWidth',1)
hold on
plot(time,chi_sh_exc1,'LineWidth',1)
hold on
plot(time,chi_pas_exc1,'m--')
grid on
title("$\chi$ Excitation 1",'Interpreter','latex')
legend("H inf","SH","Passive")
xlim([0,5])

figure
plot(time_hinf,z_hinf_exc2_1,'LineWidth',1)
hold on
plot(time,z_sh_exc2_1,'LineWidth',1)
hold on
plot(time,z_pas_exc2_1,'m--')
grid on
title("$z$ Excitation 2 (1 hz)",'Interpreter','latex')
legend("H inf","SH","Passive")
xlim([0,5])

figure
plot(time_hinf,chi_hinf_exc2_1,'LineWidth',1)
hold on
plot(time,chi_sh_exc2_1,'LineWidth',1)
hold on
plot(time,chi_pas_exc2_1,'m--')
grid on
title("$\chi$ Excitation 2 (1 hz)",'Interpreter','latex')
legend("H inf","SH","Passive")
xlim([0,5])

figure
plot(time_hinf,z_hinf_exc2_8,'LineWidth',1)
hold on
plot(time,z_sh_exc2_8,'LineWidth',1)
hold on
plot(time,z_pas_exc2_8,'m--')
grid on
title("$z$ Excitation 2 (8 hz)",'Interpreter','latex')
legend("H inf","SH","Passive")
xlim([0,5])

figure
plot(time_hinf,chi_hinf_exc2_8,'LineWidth',1)
hold on
plot(time,chi_sh_exc2_8,'LineWidth',1)
hold on
plot(time,chi_pas_exc2_8,'m--')
grid on
title("$\chi$ Excitation 2 (8 hz)",'Interpreter','latex')
legend("H inf","SH","Passive")
xlim([0,5])

%% Task 11.2
run("H_inf_ctrl_11_2_p.m")
z_hinf_exc1_t11_2_p = out11_2_p.z_exc1;
chi_hinf_exc1_t11_2_p = out11_2_p.chi_exc1;
z_hinf_exc2_1_t11_2_p = out11_2_p.z_exc2_1;
chi_hinf_exc2_1_t11_2_p = out11_2_p.chi_exc2_1;
z_hinf_exc2_8_t11_2_p = out11_2_p.z_exc2_8;
chi_hinf_exc2_8_t11_2_p = out11_2_p.chi_exc2_8;
time_hinf_11_2 = out11_2_p.tout;

figure
subplot(2,1,1)
plot(time_hinf,z_hinf_exc1,'--')
hold on
plot(time_hinf_11_2,z_hinf_exc1_t11_2_p,'LineWidth',1)
grid on
title("(11.2) $z$ Excitation 1 ",'Interpreter','latex')
legend("H inf","changed k +15%")
xlim([0,7])

subplot(2,1,2)
plot(time_hinf,chi_hinf_exc1,'--')
hold on
plot(time_hinf_11_2,chi_hinf_exc1_t11_2_p,'LineWidth',1)
grid on
title("(11.2) $\chi$ Excitation 1",'Interpreter','latex')
legend("H inf","changed k +15%")
xlim([0,7])


Fa1_1_inf = out11_2_p.Fa1_exc1;
Fa2_1_inf = out11_2_p.Fa2_exc1;
Fa1_2_1_inf = out11_2_p.Fa1_exc2_1;
Fa2_2_1_inf = out11_2_p.Fa2_exc2_1;
Fa1_2_8_inf = out11_2_p.Fa1_exc2_8;
Fa2_2_8_inf = out11_2_p.Fa2_exc2_8;

F_check_inf = [Fa1_1_inf;Fa2_1_inf;Fa1_2_1_inf;Fa2_2_1_inf;Fa1_2_8_inf;Fa2_2_8_inf];
check_magn_inf_11_2_p = max(F_check_inf)


% reduce stiffness by 15%
run("H_inf_ctrl_11_2_m.m")
z_hinf_exc1_t11_2_m = out11_2_m.z_exc1;
chi_hinf_exc1_t11_2_m = out11_2_m.chi_exc1;
z_hinf_exc2_1_t11_2_m = out11_2_m.z_exc2_1;
chi_hinf_exc2_1_t11_2_m = out11_2_m.chi_exc2_1;
z_hinf_exc2_8_t11_2_m = out11_2_m.z_exc2_8;
chi_hinf_exc2_8_t11_2_m = out11_2_m.chi_exc2_8;
time_hinf_11_2_m = out11_2_m.tout;

figure
subplot(2,1,1)
plot(time_hinf,z_hinf_exc1,'--')
hold on
plot(time_hinf_11_2_m,z_hinf_exc1_t11_2_m,'LineWidth',1)
grid on
title("(11.2) $z$ Excitation 1 ",'Interpreter','latex')
legend("H inf","changed k -15%")
xlim([0,7])

subplot(2,1,2)
plot(time_hinf,chi_hinf_exc1,'--')
hold on
plot(time_hinf_11_2_m,chi_hinf_exc1_t11_2_m,'LineWidth',1)
grid on
title("(11.2) $\chi$ Excitation 1",'Interpreter','latex')
legend("H inf","changed k -15%")
xlim([0,7])

Fa1_1_inf = out11_2_m.Fa1_exc1;
Fa2_1_inf = out11_2_m.Fa2_exc1;
Fa1_2_1_inf = out11_2_m.Fa1_exc2_1;
Fa2_2_1_inf = out11_2_m.Fa2_exc2_1;
Fa1_2_8_inf = out11_2_m.Fa1_exc2_8;
Fa2_2_8_inf = out11_2_m.Fa2_exc2_8;

F_check_inf = [Fa1_1_inf;Fa2_1_inf;Fa1_2_1_inf;Fa2_2_1_inf;Fa1_2_8_inf;Fa2_2_8_inf];
check_magn_inf_11_2_m = max(F_check_inf)

%% Task 11.3
run("H_inf_ctrl_11_3_light.m")
z_hinf_exc1_t11_3_l = out11_3_l.z_exc1;
chi_hinf_exc1_t11_3_l = out11_3_l.chi_exc1;
z_hinf_exc2_1_t11_3_l = out11_3_l.z_exc2_1;
chi_hinf_exc2_1_t11_3_l = out11_3_l.chi_exc2_1;
z_hinf_exc2_8_t11_3_l = out11_3_l.z_exc2_8;
chi_hinf_exc2_8_t11_3_l = out11_3_l.chi_exc2_8;
time_hinf_11_3_l = out11_3_l.tout;

figure
subplot(2,1,1)
plot(time_hinf,z_hinf_exc1,'--')
hold on
plot(time_hinf_11_3_l,z_hinf_exc1_t11_3_l,'LineWidth',1)
grid on
title("(11.3) $z$ Excitation 1 ",'Interpreter','latex')
legend("H inf","light car")
xlim([0,7])

subplot(2,1,2)
plot(time_hinf,chi_hinf_exc1,'--')
hold on
plot(time_hinf_11_3_l,chi_hinf_exc1_t11_3_l,'LineWidth',1')
grid on
title("(11.3) $\chi$ Excitation 1",'Interpreter','latex')
legend("H inf","light car")
xlim([0,7])


Fa1_1_inf = out11_3_l.Fa1_exc1;
Fa2_1_inf = out11_3_l.Fa2_exc1;
Fa1_2_1_inf = out11_3_l.Fa1_exc2_1;
Fa2_2_1_inf = out11_3_l.Fa2_exc2_1;
Fa1_2_8_inf = out11_3_l.Fa1_exc2_8;
Fa2_2_8_inf = out11_3_l.Fa2_exc2_8;

F_check_inf = [Fa1_1_inf;Fa2_1_inf;Fa1_2_1_inf;Fa2_2_1_inf;Fa1_2_8_inf;Fa2_2_8_inf];
check_magn_inf_11_3_light = max(F_check_inf)

run("H_inf_ctrl_11_3_heavy.m")
z_hinf_exc1_t11_3_h = out11_3_h.z_exc1;
chi_hinf_exc1_t11_3_h = out11_3_h.chi_exc1;
z_hinf_exc2_1_t11_3_h = out11_3_h.z_exc2_1;
chi_hinf_exc2_1_t11_3_h = out11_3_h.chi_exc2_1;
z_hinf_exc2_8_t11_3_h = out11_3_h.z_exc2_8;
chi_hinf_exc2_8_t11_3_h = out11_3_h.chi_exc2_8;
time_hinf_11_3_h = out11_3_h.tout;

figure
subplot(2,1,1)
plot(time_hinf,z_hinf_exc1,'--')
hold on
plot(time_hinf_11_3_h,z_hinf_exc1_t11_3_h,'LineWidth',1)
grid on
title("(11.3) $z$ Excitation 1 ",'Interpreter','latex')
legend("H inf","heavy car")
xlim([0,7])

subplot(2,1,2)
plot(time_hinf,chi_hinf_exc1,'--')
hold on
plot(time_hinf_11_3_h,chi_hinf_exc1_t11_3_h,'LineWidth',1')
grid on
title("(11.3) $\chi$ Excitation 1",'Interpreter','latex')
legend("H inf","heavy car")
xlim([0,7])

Fa1_1_inf = out11_3_h.Fa1_exc1;
Fa2_1_inf = out11_3_h.Fa2_exc1;
Fa1_2_1_inf = out11_3_h.Fa1_exc2_1;
Fa2_2_1_inf = out11_3_h.Fa2_exc2_1;
Fa1_2_8_inf = out11_3_h.Fa1_exc2_8;
Fa2_2_8_inf = out11_3_h.Fa2_exc2_8;

F_check_inf = [Fa1_1_inf;Fa2_1_inf;Fa1_2_1_inf;Fa2_2_1_inf;Fa1_2_8_inf;Fa2_2_8_inf];
check_magn_inf_11_3_heavy = max(F_check_inf)
