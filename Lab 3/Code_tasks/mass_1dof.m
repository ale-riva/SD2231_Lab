clear all
close all
clc

%% Task 1.2
mp = 0.16; %mass
cp = 0.4; %damping
kp = 6.32; %spring coeff

wn_und_1dof = sqrt(kp/mp);
c_crit = 2*sqrt(kp*mp);
r_1dof = cp/c_crit;

wn = wn_und_1dof;
run("excitations_script.m")

%% Task 1.3
s = tf('s');
Zp_Zw = (s*cp+kp)/(s^2*mp+s*cp+kp);
[peak_mag,wn_damp_1dof] = getPeakGain(Zp_Zw);


cp_und = 0;
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

%w=logspace(-2,3,2*100040); %Define frequency range
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

% dd_und = 0;
% PD_tf_und = (kp/mp)/(s^2 +dd_und/mp*s +(kp+dp)/mp);
% [peak_mag,wn_und_pd] = getPeakGain(PD_tf_und);
% wn_und_pd;
wn_und_pd = sqrt((kp+dp)/mp);


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
dd = 1.5;
dp = 0;
PD_tf = (kp/mp)/(s^2 +dd/mp*s +(kp+dp)/mp);
%r_Pd_sh = dd/(2*sqrt(kp*mp))
%To check that system is still underdamped and that wn is unchanged
[peak_mag,wn_pd] = getPeakGain( PD_tf);


%
%w=logspace(-2,3,100040); %Define frequency range
[A_PD,phi_PD]=bode(PD_tf,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A_PD(:),'LineWidth',1) %Amplitude ratio vs frequency
hold on
loglog(w,A(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')

ylim([0.01,100])
hold on
% xline(wn_und_pd,'--r','LineWidth',1)
xline(wn_und_pd,'r--','LineWidth',1)
hold on
xline(wn_und_1dof,'b--','LineWidth',1)
legend("PD Skyhook","passive damping","w_n und ","w_n und PD")

subplot(2,1,2)
semilogx(w,phi_PD(:),'LineWidth',1)
hold on 
semilogx(w,phi(:),'LineWidth',1)
grid on
hold on
% xline(wn_und_pd,'r--','LineWidth',1)
xline(wn_und_pd,'r--','LineWidth',1)
hold on
xline(wn_und_1dof,'b--','LineWidth',1)
legend("PD Skyhook","passive damping","w_n und ","w_n und PD")
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
PD_tf_mod = abs((kp/mp)./(s.^2 +dd/mp.*s +(kp+dp)/mp));
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
legend("PD ","passive damping","PD critically damped","w_n","w_n PD")
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
legend("PD ","passive damping","PD critically damped","w_n","w_n PD")
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
hi = 0; %hd*(kp+hp)/mp;
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
T = 1.5;
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
plot(time,excA,'--')
hold on
plot(time,excA_resp_pd,'LineWidth',1)
hold on
plot(time,excA_resp_pid,'LineWidth',1)
hold on
plot(time,excA_resp_sh,'LineWidth',1)
grid on
legend("input excitation","PD","PID","Skyhook")
xlim([0,7])


figure
plot(time,excB,'--')
hold on
plot(time,excB_resp_pd,'LineWidth',1)
hold on
plot(time,excB_resp_pid,'LineWidth',1)
hold on
plot(time,excB_resp_sh,'LineWidth',1)
grid on
legend("input excitation","PD","PID","Skyhook")
xlim([0,4])