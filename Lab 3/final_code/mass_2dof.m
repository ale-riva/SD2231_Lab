clear all
close all
clc

%% Task 6.1
s = tf('s');
mp = 0.16;
cp = 0.8;
kp = 6.32;
ms = 0.16;
cs = 0.05;
ks = 0.0632;

wn_und_1dof = sqrt(kp/mp);
c_crit = 2*sqrt(kp*mp);
r_1dof = cp/c_crit;

wn = wn_und_1dof;
run("excitations_script.m")

M_2dof = [ms 0; 0 mp];
C_2dof = [cs -cs;-cs (cp+cs)];
K_2dof = [ks -ks; -ks (kp+ks)];
wn_und_2dof = sqrt(eig(K_2dof,M_2dof));


%num = (cp*s+kp)*(ms*s^2+cs*s+ks);
num = (s*cp+kp)*(s*cs+ks);
den = (mp*s^2+(cp+cs)*s+(kp+ks))*(ms*s^2+cs*s+ks)-(cs*s+ks)^2
dof2_pas = num/den;

[wn_dof2,r_dof2] = damp(dof2_pas);
wn_dof2 = [wn_dof2(1);wn_dof2(3)]



%% Task 6.2
w=logspace(-2,3,2*100040); %Define frequency range

[A_dof2_pas,phi_dof2_pas]=bode(dof2_pas,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A_dof2_pas(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')
xlim([0.01,1000])
xline(wn_und_2dof,'--r','LineWidth',1)
legend("2 dof passive","wn und 2dof")

subplot(2,1,2)
semilogx(w,phi_dof2_pas(:),'LineWidth',1)
grid on
hold on
xline(wn_und_2dof,'r--','LineWidth',1)
legend("2 dof passive","wn und 2dof")
ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle("2 dof passive")


excA_resp_dof2 = lsim(dof2_pas,excA,time);
figure
plot(time,excA_resp_dof2);
grid on
hold on
plot(time,excA,'--')
legend("2 dof","excitation")
title("excA response opf 2dof passive")
xlabel("Time [s]")
ylabel("$z_p$ Amplitude",'Interpreter','latex')



excB_resp_dof2 = lsim(dof2_pas,excB,time);
figure
plot(time,excB_resp_dof2);
grid on
hold on
plot(time,excB,'--')
%ylim([-1.5e-3,3e-3])
legend("2 dof","excitation")
title("exB response opf 2dof passive")
xlabel("Time [s]")
ylabel("$z_p$ Amplitude",'Interpreter','latex')



% s = j*w_range; % complex pulsation vector
% % compute transfer function modulus 
% dof2_pas_mod = abs(((cp*s+kp).*(ms*s.^2+cs.*s+ks))./((mp*s.^2+(cp+cs)*s+(kp+ks)).*(ms*s.^2+cs*s+ks)-(cs*s+ks).^2));
% 
% 
% % PSD out 
% PSD_excC_resp_dof2 = (dof2_pas_mod.^2).*PSD_excC; 
% figure
% semilogy(w_range,PSD_excC_resp_dof2);
% grid on
% hold on
% semilogy(w_range,PSD_excC,'--')
% legend("2 dof","excitation")
% title("response to Excitation C (PSD)")
% xlabel("w range [rad/s]")
% ylabel("PSD Power Spectral Density")
% 
% 
% %step response
% y_2dof_pas = step(dof2_pas,time);
% figure
% plot(time,y_2dof_pas);
% grid on
% title("step response for 2dof passive system")
%% Task 6.3
%comparison with 1dof system
Zp_Zw = (s*cp+kp)/(s^2*mp+s*cp+kp);
[A,phi]=bode(Zp_Zw,w);


figure
subplot(2,1,1)
loglog(w,A(:),'LineWidth',1);
hold on
loglog(w,A_dof2_pas(:),'LineWidth',1);
grid on
hold on
xline(wn_und_2dof,'b--')
hold on
xline(wn_und_1dof,'r--')
legend("1 dof","2 dof","wn und 1dof","wn und 2dof");

subplot(2,1,2)
semilogx(w,phi(:),'LineWidth',1);
hold on
semilogx(w,phi_dof2_pas(:),'LineWidth',1);
grid on
hold on
xline(wn_und_2dof,'b--')
hold on
xline(wn_und_1dof,'r--')
legend("1 dof","2 dof","wn und 1dof","wn und 2dof");
sgtitle("Passive 1dof vs Passive 2dof")

figure
loglog(w,A(:),'LineWidth',1);
hold on
loglog(w,A_dof2_pas(:),'LineWidth',1);
grid on
% hold on
% xline(wn_und_1dof,'b--')
% hold on
% xline(wn_und_2dof,'r--')
legend("1 dof","2 dof");
xlim([0.1,100])
xlabel("Frequency [rad/s]")
ylabel("Magnitude")
%% Task 6.4

%FOR NATURAL FREQUENCIES: compute them with M, K by also reporting the
%equation in the report.


A_ss = [ 0 1 0 0;
      -ks/mp 0 ks/ms 0;
      0 0 0 1;
      ks/mp 0 -(kp+ks)/mp -cp/mp];
B_ss = [0 0 0;
     -1/ms 0 0;
        0 0 0;
      1/mp kp/mp cp/mp];
C_ss = [1 0 0 0];
D_ss = [0 0 0];

%inputs [F; zw ;zw_dot]
T = 0.10;
% T is underdamped --> impulse response shows oscillations. 

%just used as a check of the close loop in simulink
A_cl = [ 0 1 0 0;
      -ks/mp -T/ms ks/ms 0;
      0 0 0 1;
      ks/mp T/mp -(kp+ks)/mp -cp/mp];
B_cl = [ 0 0; 
         0 0; 
         0 0;
         kp/mp cp/mp];
C_cl = [1 0 0 0];
D_cl = [0 0];
SH_2dof_cl = ss(A_cl,B_cl,C_cl,D_cl)
SH_2dof = tf(SH_2dof_cl);

%bode of the skyhook controller with 2 dof
[A_2dof_sh,phi_2dof_sh]=bode(SH_2dof(1),w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A_2dof_sh(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
hold on
loglog(w,A_dof2_pas(:),'LineWidth',1)
ylabel('Magnitude')
% hold on
% xline(wn_dof2_sh(1),'r--','LineWidth',1)
% hold on
% xline(wn_dof2_sh(3),'o--','LineWidth',1)
% hold on
% xline(wn_und_2dof,'-.r','LineWidth',1)
% legend("sh 2dof","2 dof pas","SH wn Zs","SH wn zp","PAS wn Zs","PAs wn zp")


subplot(2,1,2)
semilogx(w,phi_2dof_sh(:),'LineWidth',1)
hold on
semilogx(w,phi_dof2_pas(:),'LineWidth',1)
grid on
% hold on
% xline(wn_dof2_sh(1),'r--','LineWidth',1)
% hold on
% xline(wn_dof2_sh(2),'o--','LineWidth',1)
% hold on
% %xline(wn_dof2(1),'-.r','LineWidth',1)
% hold on
% %xline(wn_dof2(3),'-.o','LineWidth',1)
% legend("sh 2dof","2 dof pas","SH wn Zs","SH wn zp","PAS wn Zs","PAs wn zp")
ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle("Skyhook vs Passive for 2dof")

%% Task 6.5
out2 = sim("SS_task6_5",'StartTime','0','StopTime','20','FixedStep','0.01');
time = out2.tout;
excA_resp_dof2_sh = out2.zs_2dof_excA;
excB_resp_dof2_sh = out2.zs_2dof_excB;

figure
plot(time,excA_resp_dof2_sh,'LineWidth',1);
hold on
plot(time,excA_resp_dof2,'--','LineWidth',1);
grid on
legend("2dof  sh","2dof pas")


figure
plot(time,excB_resp_dof2_sh,'LineWidth',1);
hold on
plot(time,excB_resp_dof2,'--','LineWidth',1);
grid on
legend("2dof sh","2dof pas")
%% Task 7.1 
%figures are already printed in section 6.5


