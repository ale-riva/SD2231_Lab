%% Basic excitations
%excitation a -simulink
%excition b - simulink
out  = sim('excitations','StartTime','0','StopTime','20','FixedStep','0.01');

excA = out.excA.signals.values;
excB = out.excB.signals.values;
time = out.excA.time;
% figure
% plot(time,excA)
% grid on
% hold on
% plot(time,excB)
% xlabel("Time [s]")          
% ylabel("$z_w$ Amplitude",'Interpreter','latex')
%excitation c 
syms w
S_w = (4.028*10^(-7))/(2.88*10^(-4)+0.68*w^2+w^4);
w_range = 0:0.01:25;
PSD_excC = double(subs(S_w,w,w_range));
% figure
% semilogy(w_range,PSD_excC)
% grid on
% xlabel("frequency [rad/s]")
% ylabel("Power Spectral Density")
% title("excitation 3")