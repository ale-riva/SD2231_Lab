time_trim = Time(1000:end-500)-Time(1000);

%%
% vy_mod = (SteerAngle/Ratio).*(vx_VBOX.*(lr.*(lf+lr).*Cf.*Cr-lf.*Cf.*mass.*(vx_VBOX.^2)))...
%     ./((lf+lr).^2.*Cf.*Cr+mass.*(vx_VBOX.^2).*(lr.*Cr-lf.*Cf));
% vy_mod = ((vx_VBOX(:,2).*(lr*(lf+lr)*Cf*Cr-lf*Cf*mass.*(vx_VBOX(:,2).^2)))./((lf+lr)^2*Cf*Cr+mass.*(vx_VBOX(:,2).^2)*(lr*Cr-lf*Cf))).*SteerAngle(:,2);
% 
% 
% vy_COG = vy_VBOX(:,2) + rx*yawRate_VBOX(:,2);
% yaw_acc = diff(yawRate_VBOX(time_start:time_end+1))./Ts;
% ay_COG = ay_VBOX(time_start:time_end)+rx*yaw_acc - ry*yaw_acc;
% integ = ay_VBOX(1000:end-500,2)*(1-deg2rad(rollGrad))-(yawRate_VBOX(1000:end-500,2).*vx_VBOX(1000:end-500,2));
% vy_kin = cumtrapz(time_trim,integ);
% 
% vy_COG_trim = [time_trim,vy_COG(1000:end-500)];
% vy_mod_trim = [time_trim,vy_mod(1000:end-500)];
% vy_kin_trim = [time_trim,vy_kin];
% mdl = 'washout_filter_sim';
% set_param(mdl,"StopTime",string(time_trim(end)))
% out = sim(mdl);
% vy_wf = out.vy_wf.Data;
% 
% Beta_mod = vy_mod_trim(:,2)./vx_VBOX(1000:end-500);
% Beta_kin = vy_kin_trim(:,2)./vx_VBOX(time_start:time_end);
% Beta_wf = vy_wf./vx_VBOX(time_start:time_start+size(vy_wf)-1);

%%
% figure(1)
% plot(vy_COG_trim(:,1),vy_COG_trim(:,2))
% hold on
% plot(vy_mod_trim(:,1),vy_mod_trim(:,2))
% hold on
% plot(vy_kin_trim(:,1),vy_kin_trim(:,2))
% hold on
% plot(out.vy_wf)
% xlabel("Time")
% ylabel("Speed vy")
% legend(["Ground Truth","Vy_{mod}","Vy_{kin}","Vy_{wf}"],location,'southeast')
% grid on
% title("Lateral Speed")
% 
% 
% figure(2)
% plot(time_trim,Beta_VBOX(1000:end-500))
% hold on
% plot(time_trim,Beta_mod)
% hold on
% plot(time_trim,Beta_kin)
% hold on
% plot(time_trim,Beta_wf)
% xlabel("Time")
% ylabel("Side slip amplitude")
% legend("Ground Truth","By_{mod}","By_{kin}","By_{wf}")
% title("Side slip")
% grid on

%% Used to calculate error
% disp(length(out.Betay_mod.Data))
% v1 = out.Betay_mod.Data;
% v2 = Beta_VBOX(1000:(1000+length(out.Betay_mod.Data)-1),2);
% disp(length(v1))
% disp(length(v2))


start_time_index = 1000;
beta_start = 4;
[error_mean,error_max,time_at_max,~] = errorCalc(out.Betay_mod.Data,Beta_VBOX(start_time_index:(start_time_index+length(out.Betay_mod.Data)-1),2));
disp(sprintf("Error mean for the model: %e",error_mean))
[error_mean,error_max,time_at_max,~] = errorCalc(out.Betay_kin.Data(beta_start:end),Beta_VBOX(start_time_index:(start_time_index+length(out.Betay_kin.Data(beta_start:end))-1),2));
disp(sprintf("Error mean with the kinetic: %e",error_mean))
[error_mean,error_max,time_at_max,~] = errorCalc(out.Betay_wf.Data(beta_start:end),Beta_VBOX(start_time_index:(start_time_index+length(out.Betay_kin.Data(beta_start:end))-1),2));
disp(sprintf("Error mean for the washout: %e",error_mean))
