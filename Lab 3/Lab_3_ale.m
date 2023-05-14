clear all
close all
clc



%% Task 1.2
mp = 0.16; %mass
cp = 0.4; %damping
kp = 6.32; %spring coeff

wn = sqrt(kp/mp);
r = cp/(2*sqrt(kp*mp));

%% Basic excitations
%excitation a -simulink
%excition b - simulink
out  = sim("excitations.slx");
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
Zp_Zw = (1+2*r/wn*s)/(s^2/wn^2+2*r/wn*s+1);


w=logspace(0,3,100040); %Define frequency range
[A,phi]=bode(Zp_Zw,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')
xlim([1,1000])
ylim([0.01,100])
xline(wn,'--r','LineWidth',1)
subplot(2,1,2)
semilogx(w,phi(:),'LineWidth',1)
grid on
hold on
xline(wn,'r--','LineWidth',1)
ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle("Damped System")
[peak_mag,f] = getPeakGain(Zp_Zw);
peak = mag2db(peak_mag)

%
cp_und = 0;
wn = sqrt(kp/mp);
r_undamped = cp_und/(2*sqrt(kp*mp));
Zp_Zw_undamped = (1+2*r_undamped/wn*s)/(s^2/wn^2+2*r_undamped/wn*s+1)

%w=logspace(0,3); %Define frequency range
[A_undamped,phi_undamped]=bode(Zp_Zw_undamped,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A_undamped(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')
xlim([1,1000])
ylim([0.01,1000])
xline(wn,'--r','LineWidth',1)
subplot(2,1,2)
semilogx(w,phi_undamped(:),'LineWidth',1)
grid on
hold on
xline(wn,'r--','LineWidth',1)
ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle("Undamped System")
[peak_mag,f] = getPeakGain(Zp_Zw);
peak = mag2db(peak_mag)

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
Zp_Zw_mod = abs((1+2*r/wn*s)./(s.^2/wn^2+2*r/wn*s+1));
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
r_und = r;
Zp_Zw_und = Zp_Zw;
y_und = step(Zp_Zw_und,time_s);

cp_crit = 2*sqrt(kp*mp);
r_crit = cp_crit/(2*sqrt(kp*mp));
Zp_Zw_crit = (1+2*r_crit/wn*s)/(s^2/wn^2+2*r_crit/wn*s+1);
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
PD_sh = (kp/mp)/(s^2 +dd/mp*s +(kp+dp)/mp);

w=logspace(0,3,100040); %Define frequency range
[A_PD_sh,phi_PD_sh]=bode(PD_sh,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A_PD_sh(:),'LineWidth',1) %Amplitude ratio vs frequency
hold on
loglog(w,A(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')
xlim([1,1000])
ylim([0.01,100])
xline(wn,'--r','LineWidth',1)
legend("PD Skyhook","passive damping","w_n")
subplot(2,1,2)
semilogx(w,phi_PD_sh(:),'LineWidth',1)
hold on 
semilogx(w,phi(:),'LineWidth',1)
grid on
hold on
xline(wn,'r--','LineWidth',1)
legend("PD Skyhook","passive damping","w_n")
ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle("first tuning")

%% Task 2.3

s =tf('s')
dd = 1.5*cp;
dp = 0;
PD_sh = (kp/mp)/(s^2 +dd/mp*s +(kp+dp)/mp);
r_Pd_sh = dd/(2*sqrt(kp*mp))
[wn_pd,r_pd] = damp(PD_sh);

%
w=logspace(0,3,100040); %Define frequency range
[A_PD_sh,phi_PD_sh]=bode(PD_sh,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A_PD_sh(:),'LineWidth',1) %Amplitude ratio vs frequency
hold on
loglog(w,A(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')
xlim([1,1000])
ylim([0.01,100])
xline(wn,'--r','LineWidth',1)
xline(wn_pd,'--b','LineWidth',1)
legend("PD Skyhook","passive damping","w_n","w_n PD")
subplot(2,1,2)
semilogx(w,phi_PD_sh(:),'LineWidth',1)
hold on 
semilogx(w,phi(:),'LineWidth',1)
grid on
hold on
xline(wn,'r--','LineWidth',1)
xline(wn_pd,'--b','LineWidth',1)
legend("PD Skyhook","passive damping","w_n","w_n PD")
ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle("second tuning")



excA_resp_sh = lsim(PD_sh,excA,time);
figure
plot(time,excA_resp_sh);
grid on
hold on
plot(time,excA_resp,'--')
legend("active","passive")
title("response to Excitation A (sinusoid)")
xlabel("Time [s]")
ylabel("$z_p$ Amplitude",'Interpreter','latex')


excB_resp_sh = lsim(PD_sh,excB,time);
figure
plot(time,excB_resp_sh);
grid on
hold on
plot(time,excB_resp,'--')
legend("active","passive")
title("response to Excitation B (impulse)")
xlabel("Time [s]")
ylabel("$z_p$ Amplitude",'Interpreter','latex')



s = j*w_range; % complex pulsation vector
% compute transfer function modulus 
PD_sh_mod = abs((kp/mp)./(s.^2 +dd/mp*s +(kp+dp)/mp));
% PSD out 
PSD_excC_resp_sh = (PD_sh_mod.^2).*PSD_excC; 
figure
semilogy(w_range,PSD_excC_resp_sh);
grid on
hold on
semilogy(w_range,PSD_excC_resp,'--')
legend("acrtive","passive")
title("response to Excitation C (PSD)")
xlabel("w range [rad/s]")
ylabel("PSD Power Spectral Density")

%%
%Critically damped:
s =tf('s')
dd_crit = 2*wn*mp;
PD_sh_crit = (kp/mp)/(s^2 +dd_crit/mp*s +(kp+dp)/mp);
[wn_pd_crit,r_pd_crit] = damp(PD_sh_crit);
r_Pd_sh_crit = dd_crit/(2*sqrt((kp+dp)*mp));
w=logspace(0,3,100040); %Define frequency range
[A_PD_sh_crit,phi_PD_sh_crit]=bode(PD_sh_crit,w); %Amplitude ratio (A) and phase shift (phi

w=logspace(0,3,100040); %Define frequency range
[A_PD_sh,phi_PD_sh]=bode(PD_sh,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A_PD_sh(:),'LineWidth',1) %Amplitude ratio vs frequency
hold on
loglog(w,A(:),'LineWidth',1) %Amplitude ratio vs frequency
loglog(w,A_PD_sh_crit(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')
xlim([1,1000])
ylim([0.01,100])
xline(wn,'--r','LineWidth',1)
xline(wn_pd_crit,'--b','LineWidth',1)
legend("PD Skyhook","passive damping","active critically damped","w_n","w_n PD")
subplot(2,1,2)
semilogx(w,phi_PD_sh(:),'LineWidth',1)
hold on 
semilogx(w,phi(:),'LineWidth',1)
hold on
semilogx(w,phi_PD_sh_crit(:),'LineWidth',1)
hold on
xline(wn,'r--','LineWidth',1)
xline(wn_pd_crit,'--b','LineWidth',1)
legend("PD Skyhook","passive damping","active critically damped","w_n","w_n PD")
grid on

ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle("skyhook critically damped")

% As a double check, the step response of the system should be with no
% oscillations and the fastest to converge
y_pd_sh_crit = step(PD_sh_crit,time);
dd_t1 = 0.9*2*wn*mp;
dd_t2 = 1.1*2*wn*mp;
PD_sh_t1 = (kp/mp)/(s^2 +dd_t1/mp*s +(kp+dp)/mp);
PD_sh_t2 = (kp/mp)/(s^2 +dd_t2/mp*s +(kp+dp)/mp);
PD_sh_crit = (kp/mp)/(s^2 +dd_crit/mp*s +(kp+dp)/mp);


y_pd_sh_t1 = step(PD_sh_t1,time)
y_pd_sh_t2 = step(PD_sh_t2,time)

figure
plot(time,y_pd_sh_crit)
hold on
plot(time,y_pd_sh_t1)
hold on
plot(time,y_pd_sh_t2)
grid on
legend("$d_d = c_{critical}$","$d_d = 0.9c_{critical}$","$d_d = 1.1c_{critical}$",'Interpreter','latex')
ylim([0.99,1.01])
xlim([0,4])
title("check for critically damped active system")
%% Task 3.1
close all
s = tf('s')
hd = 1.5*cp;
hp = 0;
hi = 100;
%wn_pid = (kp+;
PID_sh = (kp)/(s^2*mp + s*hd + (kp+hp)+hi/s);
[wn_pid,r_pid] = damp(PID_sh)
w=logspace(0,3,100040); %Define frequency range
[A_PID_sh,phi_PID_sh]=bode(PID_sh,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A_PID_sh(:),'LineWidth',1) %Amplitude ratio vs frequency
hold on
loglog(w,A(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')
xlim([1,1000])
ylim([0.01,100])
xline(wn,'--r','LineWidth',1)
xline(wn_pid(end),'b--','LineWidth',1)
legend("PID Skyhook","passive damping","w_n",'w_n PID')
subplot(2,1,2)
semilogx(w,phi_PID_sh(:),'LineWidth',1)
hold on 
semilogx(w,phi(:),'LineWidth',1)
grid on
hold on
xline(wn,'r--','LineWidth',1)
xline(wn_pid(end),'b--','LineWidth',1)
legend("PID Skyhook","passive damping","w_n",'w_n PID')
ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle("")

%comparison with excitations
excA_resp_sh = lsim(PID_sh,excA,time);
figure
plot(time,excA_resp_sh);
grid on
hold on
plot(time,excA_resp,'--')
legend("PID","passive")
title("response to Excitation A (sinusoid)")
xlabel("Time [s]")
ylabel("$z_p$ Amplitude",'Interpreter','latex')


excB_resp_sh = lsim(PID_sh,excB,time);
figure
plot(time,excB_resp_sh);
grid on
hold on
plot(time,excB_resp,'--')
legend("PID","passive")
title("response to Excitation B (impulse)")
xlabel("Time [s]")
ylabel("$z_p$ Amplitude",'Interpreter','latex')



s = j*w_range; % complex pulsation vector
% compute transfer function modulus 
PID_sh_mod = abs((kp*s)./(s.^3*mp + s.^2*hd + (kp+hp)*s+hi));
% PSD out 
PSD_excC_resp_sh = (PID_sh_mod.^2).*PSD_excC; 
figure
semilogy(w_range,PSD_excC_resp_sh);
grid on
hold on
semilogy(w_range,PSD_excC_resp,'--')
legend("PID","passive")
title("response to Excitation C (PSD)")
xlabel("w range [rad/s]")
ylabel("PSD Power Spectral Density")


