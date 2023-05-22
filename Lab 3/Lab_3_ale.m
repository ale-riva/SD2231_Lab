clear all
close all
clc

%% Task 1.2
mp = 0.16; %mass
cp = 0.4; %damping
kp = 6.32; %spring coeff

wn_und_1dof = sqrt(kp/mp);
r_1dof = cp/(2*sqrt(kp*mp));

%% Basic excitations
%excitation a -simulink
%excition b - simulink
out  = sim('excitations','StartTime','0','StopTime','20','FixedStep','0.01');

excA = out.excA.signals.values;
excB = out.excB.signals.values;
time = out.excA.time;
figure
plot(time,excA)
grid on
hold on
plot(time,excB)
xlabel("Time [s]")          
ylabel("$z_w$ Amplitude",'Interpreter','latex')
%excitation c 
syms w
S_w = (4.028*10^(-7))/(2.88*10^(-4)+0.68*w^2+w^4);
w_range = 0:0.01:25;
PSD_excC = double(subs(S_w,w,w_range));
figure
semilogy(w_range,PSD_excC)
grid on
xlabel("frequency [rad/s]")
ylabel("Power Spectral Density")
title("excitation 3")


%% Task 1.3
s = tf('s');
Zp_Zw = (s*cp+kp)/(s^2*mp+s*cp+kp);
[peak_mag,wn_damp_1dof] = getPeakGain(Zp_Zw);

cp_und = 0;
wn = sqrt(kp/mp);
r_undamped = cp_und/(2*sqrt(kp*mp));
Zp_Zw_undamped = (1+2*r_undamped/wn*s)/(s^2/wn^2+2*r_undamped/wn*s+1)
w=logspace(-2,3,2*100040); %Define frequency range
[A_undamped,phi_undamped]=bode(Zp_Zw_undamped,w);


[A,phi]=bode(Zp_Zw,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
xlim([0.01,1000])
loglog(w,A(:),'r','LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')

xline(wn_damp_1dof,'--r','LineWidth',1)
legend("damped","wn damped")

subplot(2,1,2)
semilogx(w,phi(:),'r','LineWidth',1)
hold on
xline(wn_damp_1dof(1),'r--','LineWidth',1)
grid on
legend("damped","wn damped")
ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle("Damped System")

w=logspace(-2,3,2*100040); %Define frequency range
[A,phi]=bode(Zp_Zw,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
xlim([0.01,1000])
loglog(w,A(:),'r','LineWidth',1) %Amplitude ratio vs frequency
hold on
loglog(w,A_undamped(:),'b','LineWidth',1)
grid on
ylabel('Magnitude')
xlim([1,70])
xline(wn_und_1dof,'--b','LineWidth',1)
hold on
xline(wn_damp_1dof,'--r','LineWidth',1)
legend("damped","undamped","wn undamped","wn damped")

subplot(2,1,2)
semilogx(w,phi(:),'r','LineWidth',1)
hold on
semilogx(w,phi_undamped(:),'b','LineWidth',1)
grid on
hold on
xline(wn_und_1dof,'b--','LineWidth',1)
hold on
xline(wn_damp_1dof(1),'r--','LineWidth',1)
xlim([1,70])
legend("damped","undamped","wn undamped","wn damped")
ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle("Damped System vs Undamped System")

%% Task 1.4
excA_resp = lsim(Zp_Zw,excA,time);
figure
plot(time,excA_resp);
grid on
hold on
plot(time,excA,'--')
legend("response","excitation")
title("response to Excitation A (sinusoid)")
xlabel("Time [s]")
ylabel("$z_p$ Amplitude",'Interpreter','latex')


excB_resp = lsim(Zp_Zw,excB,time);
figure
plot(time,excB_resp);
grid on
hold on
plot(time,excB,'--')
legend("response","excitation")
title("response to Excitation B (impulse)")
xlabel("Time [s]")
ylabel("$z_p$ Amplitude",'Interpreter','latex')

s = j*w_range; % complex pulsation vector
% compute transfer function modulus 
Zp_Zw_mod = abs((1+2*r_1dof/wn_und_1dof*s)./(s.^2/wn_und_1dof^2+2*r_1dof/wn_und_1dof*s+1));
% PSD out 
PSD_excC_resp = (Zp_Zw_mod.^2).*PSD_excC; 
figure
semilogy(w_range,PSD_excC_resp);
grid on
hold on
semilogy(w_range,PSD_excC,'--')
legend("response","excitation")
title("response to Excitation C (PSD)")
xlabel("w range [rad/s]")
ylabel("PSD Power Spectral Density")

%Step response for underdamped, critically damped, overdamped syst.
s= tf('s');

time_s = [0:0.01:6];
cp_und = cp;
r_und = r_1dof;
Zp_Zw_und = Zp_Zw;
y_und = step(Zp_Zw_und,time_s);

cp_crit = 2*sqrt(kp*mp);
r_crit = cp_crit/(2*sqrt(kp*mp));
Zp_Zw_crit = (1+2*r_crit/wn_und_1dof*s)/(s^2/wn_und_1dof^2+2*r_crit/wn_und_1dof*s+1);
y_crit = step(Zp_Zw_crit,time_s);

cp_over = 1.5*2*sqrt(kp*mp);
r_over = cp_over/(2*sqrt(kp*mp));
Zp_Zw_over = (1+2*r_over/wn*s)/(s^2/wn^2+2*r_over/wn*s+1);
y_over = step(Zp_Zw_over,time_s);

figure
plot(time_s,y_und)
hold on
plot(time_s,y_crit)
hold on
plot(time_s,y_over)
legend("$\zeta <1$","$\zeta = 1$","$\zeta >1$",'Interpreter','latex')
grid on
title("step response for different $\zeta$",'Interpreter','latex')

%% Task 2.1
dd = cp;
dp = 0;
 PD_tf = (kp/mp)/(s^2 +dd/mp*s +(kp+dp)/mp);
[peak_mag,wn_pd] = getPeakGain( PD_tf);

dd_und = 0;
PD_tf_und = (kp/mp)/(s^2 +dd_und/mp*s +(kp+dp)/mp);
[peak_mag,wn_und_pd] = getPeakGain(PD_tf_und);
wn_und_pd;



%w=logspace(0,3,100040); %Define frequency range
[A_PD,phi_PD]=bode(PD_tf,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A_PD(:),'LineWidth',1) %Amplitude ratio vs frequency
hold on
loglog(w,A(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')
ylim([0.01,100])
xline(wn_und_pd,'--r','LineWidth',1)
hold on
xline(wn_pd,'--b','LineWidth',1)

legend("PD ","passive damping","w_n und","w_n pd")

subplot(2,1,2)
semilogx(w,phi_PD(:),'LineWidth',1)
hold on 
semilogx(w,phi(:),'LineWidth',1)
grid on
hold on
xline(wn_und_pd,'r--','LineWidth',1)
hold on
xline(wn_pd,'--b','LineWidth',1)
legend("PD ","passive damping","w_n und","w_n pd")
ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle("first tuning")

%% Task 2.3 PD

s =tf('s')
dd = 2*cp;
dp = 0;
PD_tf = (kp/mp)/(s^2 +dd/mp*s +(kp+dp)/mp);
%r_Pd_sh = dd/(2*sqrt(kp*mp))
%To check that system is still underdamped and that wn is unchanged
%[wn_pd,r_pd] = damp(PD_tf);


%
%w=logspace(0,3,100040); %Define frequency range
[A_PD,phi_PD]=bode(PD_tf,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A_PD(:),'LineWidth',1) %Amplitude ratio vs frequency
hold on
loglog(w,A(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')

ylim([0.01,100])
% xline(wn_und_pd,'--r','LineWidth',1)
xline(wn_pd,'--b','LineWidth',1)
legend("PD Skyhook","passive damping","w_n","w_n PD")

subplot(2,1,2)
semilogx(w,phi_PD(:),'LineWidth',1)
hold on 
semilogx(w,phi(:),'LineWidth',1)
grid on
hold on
% xline(wn_und_pd,'r--','LineWidth',1)
xline(wn_pd,'--b','LineWidth',1)
legend("PD Skyhook","passive damping","w_n","w_n PD")
ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle("second tuning")

excA_resp_pd = lsim(PD_tf,excA,time);
figure
plot(time,excA_resp_pd);
grid on
hold on
plot(time,excA_resp,'--')
legend("active","passive")
title("response to Excitation A (sinusoid)")
xlabel("Time [s]")
ylabel("$z_p$ Amplitude",'Interpreter','latex')


excB_resp_pd = lsim(PD_tf,excB,time);
figure
plot(time,excB_resp_pd);
grid on
hold on
plot(time,excB_resp,'--')
legend("active","passive")
title("response to Excitation B (impulse)")
xlabel("Time [s]")
ylabel("$z_p$ Amplitude",'Interpreter','latex')



s = j*w_range; % complex pulsation vector
% compute transfer function modulus 
PD_tf_mod = abs((kp/mp)./(s.^2 +dd/mp*s +(kp+dp)/mp));
% PSD out 
PSD_excC_resp_pd = (PD_tf_mod.^2).*PSD_excC; 
figure
semilogy(w_range,PSD_excC_resp_pd);
grid on
hold on
semilogy(w_range,PSD_excC_resp,'--')
legend("acrtive","passive")
title("response to Excitation C (PSD)")
xlabel("w range [rad/s]")
ylabel("PSD Power Spectral Density")

%% Extra task 2
%Critically damped:
s =tf('s')
dd_crit = 2*wn_und_pd*mp;
PD_tf_crit = (kp/mp)/(s^2 +dd_crit/mp*s +(kp+dp)/mp);
%[wn_pd_crit,r_pd_crit] = damp(PD_tf_crit);
r_Pd_sh_crit = dd_crit/(2*sqrt((kp+dp)*mp));
%w=logspace(0,3,100040); %Define frequency range
[A_PD_crit,phi_PD_crit]=bode(PD_tf_crit,w); %Amplitude ratio (A) and phase shift (phi

%w=logspace(0,3,100040); %Define frequency range
[A_PD,phi_PD]=bode(PD_tf,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A_PD(:),'LineWidth',1) %Amplitude ratio vs frequency
hold on
loglog(w,A(:),'LineWidth',1) %Amplitude ratio vs frequency
loglog(w,A_PD_crit(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')

ylim([0.01,100])
xline(wn_und_1dof,'--r','LineWidth',1)
% xline(wn_und_pd,'--b','LineWidth',1)
xline(wn_pd,'--b','LineWidth',1)
legend("PD Skyhook","passive damping","active critically damped","w_n","w_n PD")
subplot(2,1,2)
semilogx(w,phi_PD(:),'LineWidth',1)
hold on 
semilogx(w,phi(:),'LineWidth',1)
hold on
semilogx(w,phi_PD_crit(:),'LineWidth',1)
hold on
xline(wn_und_1dof,'r--','LineWidth',1)
% xline(wn_und_pd,'--b','LineWidth',1)
xline(wn_pd,'--b','LineWidth',1)
legend("PD Skyhook","passive damping","active critically damped","w_n","w_n PD")
grid on

ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle("PD critically damped")

% As a double check, the step response of the system should be with no
% oscillations and the fastest to converge
y_pd_sh_crit = step(PD_tf_crit,time);
dd_t1 = 0.9*2*wn*mp;
dd_t2 = 1.1*2*wn*mp;
PD_tf_t1 = (kp/mp)/(s^2 +dd_t1/mp*s +(kp+dp)/mp);
PD_tf_t2 = (kp/mp)/(s^2 +dd_t2/mp*s +(kp+dp)/mp);
PD_tf_crit = (kp/mp)/(s^2 +dd_crit/mp*s +(kp+dp)/mp);


y_pd_t1 = step(PD_tf_t1,time)
y_pd_t2 = step(PD_tf_t2,time)

figure
plot(time,y_pd_sh_crit)
hold on
plot(time,y_pd_t1)
hold on
plot(time,y_pd_t2)
grid on
legend("$d_d = c_{critical}$","$d_d = 0.9c_{critical}$","$d_d = 1.1c_{critical}$",'Interpreter','latex')
ylim([0.99,1.01])
xlim([0,4])
title("check for critically damped active system")
%% Task 3.1 PID

s = tf('s');
hd = 1.5;%cp;
hp = 0;
hi = 0.2; %hd*(kp+hp)/mp;
PID_tf = (kp*s)/(s^3*mp + s^2*hd + (kp+hp)*s+hi);

hd_und = 0;
PID_und_tf = (kp*s)/(s^3*mp + s^2*hd_und + (kp+hp)*s+hi);
[peak_mag,wn_und_pid] = getPeakGain(PID_und_tf);



[A_PID,phi_PID]=bode(PID_tf,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A_PID(:),'LineWidth',1) %Amplitude ratio vs frequency
hold on
loglog(w,A(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')
xlim([0.01,1000])
ylim([0.01,100])
xline(wn_und_1dof,'--r','LineWidth',1)
hold on
xline(wn_und_pid,'b--','LineWidth',1)
legend("PID Skyhook","passive damping","w_n und",'w_n und PID')

subplot(2,1,2)
semilogx(w,phi_PID(:),'LineWidth',1)
hold on 
semilogx(w,phi(:),'LineWidth',1)
grid on
hold on
xline(wn_und_1dof,'--r','LineWidth',1)
hold on
xline(wn_und_pid,'b--','LineWidth',1)
legend("PID Skyhook","passive damping","w_n und",'w_n und PID')
ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle(hi)

excA_resp_pid = lsim(PID_tf,excA,time);
figure
plot(time,excA_resp_pid);
grid on
hold on
plot(time,excA_resp,'--')
legend("PID","passive")
title("response to Excitation A (sinusoid)")
xlabel("Time [s]")
ylabel("$z_p$ Amplitude",'Interpreter','latex')


excB_resp_pid = lsim(PID_tf,excB,time);
figure
plot(time,excB_resp_pid);
grid on
hold on
plot(time,excB_resp,'--')
legend("PID","passive")
title("response to Excitation B (impulse)")
xlabel("Time [s]")
ylabel("$z_p$ Amplitude",'Interpreter','latex')


s = j*w_range; % complex pulsation vector
% compute transfer function modulus 
PID_tf_mod = abs((kp*s)./(s.^3*mp + s.^2*hd + (kp+hp)*s+hi));
% PSD out 
PSD_excC_resp_pid = (PID_tf_mod.^2).*PSD_excC; 
figure
semilogy(w_range,PSD_excC_resp_pid);
grid on
hold on
semilogy(w_range,PSD_excC_resp,'--')
legend("PID","passive")
title("response to Excitation C (PSD)")
xlabel("w range [rad/s]")
ylabel("PSD Power Spectral Density")

%step response
y_pid = step(PID_tf,time);
figure
plot(time,y_pid)
grid on
title("step response for PID controlled system")


%% Task 4.1 Skyhook control
s = tf('s')
T = 3*cp;
SH_ctrl = kp/(mp*s^2+T*s+kp);
%[wn_sh,r_sh] = damp(SH_ctrl);
SH_und_ctrl = kp/(mp*s^2+kp);
[peak_mag,wn_und_sh] = getPeakGain(SH_und_ctrl);



%w=logspace(0,3,100040); %Define frequency range
[A_sh,phi_sh]=bode(SH_ctrl,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A_sh(:),'LineWidth',1) %Amplitude ratio vs frequency
hold on
loglog(w,A(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')

ylim([0.01,100])
xline(wn,'--r','LineWidth',1)
xline(wn_und_sh,'b--','LineWidth',1)
legend("Skyhook ","passive damping","w_n",'w_n sh')
subplot(2,1,2)
semilogx(w,phi_sh(:),'LineWidth',1)
hold on 
semilogx(w,phi(:),'LineWidth',1)
grid on
hold on
xline(wn,'r--','LineWidth',1)
xline(wn_und_sh,'b--','LineWidth',1)
legend("Skyhook","passive damping","w_n",'w_n sh')
ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')

excA_resp_sh = lsim(SH_ctrl,excA,time);
figure
plot(time,excA_resp_sh);
grid on
hold on
plot(time,excA_resp,'--')
legend("Skyhook","passive")
title("response to Excitation A (sinusoid)")
xlabel("Time [s]")
ylabel("$z_p$ Amplitude",'Interpreter','latex')


excB_resp_sh = lsim(SH_ctrl,excB,time);
figure
plot(time,excB_resp_sh);
grid on
hold on
plot(time,excB_resp,'--')
legend("Skyhook","passive")
title("response to Excitation B (impulse)")
xlabel("Time [s]")
ylabel("$z_p$ Amplitude",'Interpreter','latex')

s = j*w_range; % complex pulsation vector
% compute transfer function modulus 
SH_mod = abs(kp./(mp*s.^2+T*s+kp));
% PSD out 
PSD_excC_resp_sh = (SH_mod.^2).*PSD_excC; 
figure
semilogy(w_range,PSD_excC_resp_sh);
grid on
hold on
semilogy(w_range,PSD_excC_resp,'--')
legend("Skyhook","passive")
title("response to Excitation C (PSD)")
xlabel("w range [rad/s]")
ylabel("PSD Power Spectral Density")

%step response
y_sh = step(SH_ctrl,time);
figure
plot(time,y_sh)
grid on
title("step response for Skyhook controlled system")



%% Task 5.1
figure
loglog(w,A_undamped(:), 'LineWidth',1)
hold on
loglog(w,A(:),'LineWidth',1) %passive damped
hold on
loglog(w,A_PD(:),'LineWidth',1); %PD
hold on
loglog(w,A_PID(:),'LineWidth',1); %PID
hold on
loglog(w,A_sh(:),'LineWidth',1); % skyhook
legend("undamped","passive","PD","PID","Skyhook")
grid on
xlabel("frequency [rad/s]")
ylim([0.01,1e5])
ylabel("amplitude")

%% Task 5.2
figure
plot(time,excA_resp_pd)
hold on
plot(time,excA_resp_pid)
hold on
plot(time,excA_resp_sh)
grid on
legend("PD","PID","Skyhook")
xlim([0,7])


figure
plot(time,excB_resp_pd)
hold on
plot(time,excB_resp_pid)
hold on
plot(time,excB_resp_sh)
grid on
legend("PD","PID","Skyhook")
xlim([0,4])
%% Task 6.1
s = tf('s');
mp = 0.16;
cp = 0.8;
kp = 6.32;
ms = 0.16;
cs = 0.05;
ks = 0.0632;

%num = (cp*s+kp)*(ms*s^2+cs*s+ks);
num = (s*cp+kp)*(s*cs+ks);
den = (mp*s^2+(cp+cs)*s+(kp+ks))*(ms*s^2+cs*s+ks)-(cs*s+ks)^2
dof2_pas = num/den;

% [wn_dof2,r_dof2] = damp(dof2_pas);
% wn_dof2 = [wn_dof2(1);wn_dof2(3)];

cp_und = 0;cs_und = 0;
num_und = (s*cp_und+kp)*(s*cs_und+ks);
den_und = (mp*s^2+(cp_und+cs_und)*s+(kp+ks))*(ms*s^2+cs_und*s+ks)-(cs_und*s+ks)^2
dof2_und = num_und/den_und;
[wn_und_dof2,r_und_dof2] = damp(dof2_und)
wn_und_dof2 = [wn_und_dof2(1);wn_und_dof2(3)];


w=logspace(-2,3,2*100040); %Define frequency range
[A_dof2_und,phi_dof2_und]=bode(dof2_und,w); %Amplitude ratio (A) and phase shift (phi)

[A_dof2_pas,phi_dof2_pas]=bode(dof2_pas,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A_dof2_pas(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')
xlim([0.01,1000])
xline(wn_und_dof2,'--r','LineWidth',1)
legend("2 dof passive","wn und 2dof")

subplot(2,1,2)
semilogx(w,phi_dof2_pas(:),'LineWidth',1)
grid on
hold on
xline(wn_und_dof2,'r--','LineWidth',1)
legend("2 dof passive","wn und 2dof")
ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle("2 dof passive")

%% 
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



s = j*w_range; % complex pulsation vector
% compute transfer function modulus 
dof2_pas_mod = abs(((cp*s+kp).*(ms*s.^2+cs.*s+ks))./((mp*s.^2+(cp+cs)*s+(kp+ks)).*(ms*s.^2+cs*s+ks)-(cs*s+ks).^2));


% PSD out 
PSD_excC_resp_dof2 = (dof2_pas_mod.^2).*PSD_excC; 
figure
semilogy(w_range,PSD_excC_resp_dof2);
grid on
hold on
semilogy(w_range,PSD_excC,'--')
legend("2 dof","excitation")
title("response to Excitation C (PSD)")
xlabel("w range [rad/s]")
ylabel("PSD Power Spectral Density")


%step response
y_sh = step(dof2_pas,time)
figure
plot(time,y_sh)
grid on
title("step response for 2dof passive system")
%% Task 6.3
figure
subplot(2,1,1)
loglog(w,A(:),'LineWidth',1);
hold on
loglog(w,A_dof2_pas(:),'LineWidth',1);
grid on
hold on
xline(wn_und_dof2,'b--')
hold on
xline(wn_und_1dof,'r--')
legend("1 dof","2 dof","wn 1dof","wn 2dof");

subplot(2,1,2)
semilogx(w,phi(:),'LineWidth',1);
hold on
semilogx(w,phi_dof2_pas(:),'LineWidth',1);
grid on
hold on
xline(wn_und_1dof,'b--')
hold on
xline(wn_dof2,'r--')
legend("1 dof","2 dof","wn 1dof","wn 2dof");
sgtitle("Passive 1dof vs Passive 2dof")

figure
loglog(w,A(:),'LineWidth',1);
hold on
loglog(w,A_dof2_pas(:),'LineWidth',1);
grid on
hold on
xline(wn_und_1dof,'b--')
hold on
xline(wn_und_dof2,'r--')
legend("1 dof","2 dof","wn 1dof","wn 2dof");
xlim([0.1,100])
xlabel("Frequency [rad/s]")
%% Task 6.4

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
T = 1.5;
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
[wn_dof2_sh,r_dof2_sh] = damp(SH_2dof_cl);


%bode of the skyhook controller with 2 dof
[A_2dof_sh,phi_2dof_sh]=bode(SH_2dof(1),w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A_2dof_sh(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
hold on
loglog(w,A_dof2_pas(:),'LineWidth',1)
ylabel('Magnitude')
hold on
xline(wn_dof2_sh(1),'r--','LineWidth',1)
hold on
xline(wn_dof2_sh(3),'o--','LineWidth',1)
hold on
xline(wn_dof2(1),'-.r','LineWidth',1)
hold on
xline(wn_dof2(3),'-.o','LineWidth',1)
legend("sh 2dof","2 dof pas","SH wn Zs","SH wn zp","PAS wn Zs","PAs wn zp")


subplot(2,1,2)
semilogx(w,phi_2dof_sh(:),'LineWidth',1)
hold on
semilogx(w,phi_dof2_pas(:),'LineWidth',1)
grid on
hold on
xline(wn_dof2_sh(1),'r--','LineWidth',1)
hold on
xline(wn_dof2_sh(3),'o--','LineWidth',1)
hold on
xline(wn_dof2(1),'-.r','LineWidth',1)
hold on
xline(wn_dof2(3),'-.o','LineWidth',1)
legend("sh 2dof","2 dof pas","SH wn Zs","SH wn zp","PAS wn Zs","PAs wn zp")
ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle("Skyhook vs Passive for 2dof")


out2 = sim("SS_task6_5",'StartTime','0','StopTime','20','FixedStep','0.01');
time = out2.tout;
excA_resp_dof2_sh = out2.zs_2dof_excA;
excB_resp_dof2_sh = out2.zs_2dof_excB;

figure
plot(time,excA_resp_dof2_sh);
hold on
plot(time,excA_resp_dof2);
grid on
legend("2dof  sh","2dof pas")
title("excA resp Skyhook vs Passive 2dof")

figure
plot(time,excB_resp_dof2_sh);
hold on
plot(time,excB_resp_dof2);
grid on
legend("2dof sh","2dof pas")
title("excB resp Skyhook vs Passive 2dof")


%% Task 8.1 
m = 22000; J = 700000;
c = 40000;
k = 600000;
k1 = k; k2= k;
syms c1 c2
L = 6;  L1= L; L2= L;

A_ss_2dof_vh_sym = [0  1 0 0;
    -(k1+k2)/m ,-(c1+c2)/m,(k1*L1-k2*L2)/m,(c1*L1-c2*L2)/m;
     0, 0, 0 ,1;
     (k1*L1-k2*L2)/J, (c1*L1-c2*L2)/J,(-k1*L1^2-k2*L2^2)/J,(-c1*L1^2-c2*L2^2)/J]
B_ss_2dof_vh_sym = [0 0 0 0;
                k1/m c1/m  k2/m  c2/m;
                 0     0     0     0;
                 -(k1*L1)/m  -c1*L1/J k2*L2/J c2*L2/J]
C_ss_2dof_vh = [1 0 0 0;
                0 0 1 0];
D_ss_2dof_vh = [0 0 0 0; 0 0 0 0];

A_ss_2dof_vh = eval(subs(A_ss_2dof_vh_sym,{c1,c2},[c,c]));
B_ss_2dof_vh = eval(subs(B_ss_2dof_vh_sym,{c1,c2},[c,c]));
out8_2 = sim("SS_task8_2",'StartTime','0','StopTime','20','FixedStep','0.01');

time = out8_2.tout;
exc1_z_resp_2dof_vh = out8_2.z_exc1;
exc1_chi_resp_2dof_vh = out8_2.chi_exc1;
exc2_8_z_resp_2dof_vh = out8_2.z_exc2_8;
exc2_8_chi_resp_2dof_vh = out8_2.chi_exc2_8;
exc2_1_z_resp_2dof_vh = out8_2.z_exc2_1;
exc2_1_chi_resp_2dof_vh = out8_2.chi_exc2_1;
%% Task 8.2
A_ss_2dof_vh_und = eval(subs(A_ss_2dof_vh_sym,{c1,c2},[0,0]));
B_ss_2dof_vh_und = eval(subs(B_ss_2dof_vh_sym,{c1,c2},[0,0]));
veh_2dof_tf_und = tf(ss(A_ss_2dof_vh_und,B_ss_2dof_vh_und,C_ss_2dof_vh,D_ss_2dof_vh));
veh_2dof_tf_und.InputName = ["z_{w1}"," z_{w1} dot","z_{w2}"," z_{w2} dot"];
veh_2dof_tf_und.OutputName = ["zeta","chi"];
%%
veh_2dof_tf = tf(ss(A_ss_2dof_vh,B_ss_2dof_vh,C_ss_2dof_vh,D_ss_2dof_vh));
veh_2dof_tf.InputName = ["z_{w1}"," z_{w1} dot","z_{w2}"," z_{w2} dot"];
veh_2dof_tf.OutputName = ["zeta","chi"];
[A_veh_2dof,phi_veh_2dof]=bode(veh_2dof_tf,w); 
[A_veh_2dof_und,phi_veh_2dof_und]=bode(veh_2dof_tf_und,w); 
%%
[wn_veh,r_veh] = damp(veh_2dof_tf)
mag_lin_z_zw1 = squeeze(A_veh_2dof_und(1,1,:)); % Convert magnitude to linear scale
mag_lin_chi_zw1 = squeeze(A_veh_2dof_und(2,1,:));
[~, peak_indices_z_zw1] = findpeaks(mag_lin_z_zw1); % Find peak indices
[~, peak_indices_chi_zw1] = findpeaks(mag_lin_chi_zw1); % Find peak indices
peak_indices=[peak_indices_z_zw1,peak_indices_chi_zw1]
nat_freq = w(peak_indices) % Natural frequencies



tmp1 = reshape(A_veh_2dof_und(1,1,:),[],size(A_veh_2dof_und(1,1,:), 3));
tmp2 = reshape(A_veh_2dof_und(2,1,:),[],size(A_veh_2dof_und(2,1,:), 3));

figure
loglog(w,tmp1)
hold on
loglog(w,tmp2)
grid on

hold on
xline(nat_freq(1))
hold on
xline(nat_freq(2))
legend("z zw1","chi zw1","wn z","wn chi");
xlim([5,15])
%%
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

