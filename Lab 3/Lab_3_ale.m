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

time = [0:0.01:6];
cp_und = cp;
r_und = r;
Zp_Zw_und = Zp_Zw;
y_und = step(Zp_Zw_und,time);

cp_crit = 2*sqrt(kp*mp);
r_crit = cp_crit/(2*sqrt(kp*mp));
Zp_Zw_crit = (1+2*r_crit/wn*s)/(s^2/wn^2+2*r_crit/wn*s+1);
y_crit = step(Zp_Zw_crit,time);

cp_over = 1.5*2*sqrt(kp*mp);
r_over = cp_over/(2*sqrt(kp*mp));
Zp_Zw_over = (1+2*r_over/wn*s)/(s^2/wn^2+2*r_over/wn*s+1);
y_over = step(Zp_Zw_over,time);

figure
plot(time,y_und)
hold on
plot(time,y_crit)
hold on
plot(time,y_over)
legend("$\zeta <1$","$\zeta = 1$","$\zeta >1$",'Interpreter','latex')
grid on
title("step response for different $\zeta$",'Interpreter','latex')

