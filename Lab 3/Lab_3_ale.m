clear all
close all
clc

%% Task 1.2
mp = 0.16; %mass
cp = 0.4; %damping
kp = 6.32; %spring coeff

wn = sqrt(kp/mp);
r = cp/(2*sqrt(kp*mp));

%% Task 1.3
s = tf('s');
Zp_Zw = (1+2*r/wn*s)/(s^2/wn^2+2*r/wn*s+1);


w=logspace(0,3); %Define frequency range
[A,phi]=bode(Zp_Zw,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')
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
cp = 0;
wn = sqrt(kp/mp);
r = cp/(2*sqrt(kp*mp));
Zp_Zw_undamped = (1+2*r/wn*s)/(s^2/wn^2+2*r/wn*s+1)

w=logspace(0,3); %Define frequency range
[A,phi]=bode(Zp_Zw_undamped,w); %Amplitude ratio (A) and phase shift (phi)
figure
subplot(2,1,1)
loglog(w,A(:),'LineWidth',1) %Amplitude ratio vs frequency
grid on
ylabel('Magnitude')
ylim([0.01,100])
xline(wn,'--r','LineWidth',1)
subplot(2,1,2)
semilogx(w,phi(:),'LineWidth',1)
grid on
hold on
xline(wn,'r--','LineWidth',1)
ylabel('Phase angle(Degrees)')
xlabel('Frequency(rad/s)')
sgtitle("Undamped system")
