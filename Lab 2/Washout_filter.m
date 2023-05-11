time_trim = Time(1000:end-500)-Time(1000);



Beta_VBOX_smooth=smooth(Beta_VBOX(:,2),0.01,'rlowess');



% Portion where to evaluate errors of the estimation (+time_start)
start_time_index = 1001;
time_start = Time(start_time_index);


sim1 = [4 37.28];
sim2 = [5.5 21.78];
sim3 = [7 11.510];
sim4 = [0 27.23];

switch vbox_file_name
    case 'S90__036.VBO'
        sim_txt = "sim1";
        sim_ok = sim1;
    case 'S90__038.VBO'
        sim_txt = "sim2";
        sim_ok = sim2;
    case 'S90__040.VBO'
        sim_txt = "sim3";
        sim_ok = sim3;
    case 'S90__041.VBO'
        sim_txt = "sim4";
        sim_ok = sim4;
end
sim_ok
sim = sim_ok.*100+1;
sim_truth = (sim_ok +time_start).*100+1;

figure
plot(out.Betay_mod.Data(sim(1):sim(2)))
hold on
plot(out.Betay_kin.Data(sim(1):sim(2)))
hold on
plot(out.Betay_wf.Data(sim(1):sim(2)))
hold on
plot(out.Betay_wf_var.Data(sim(1):sim(2)))
hold on
plot(Beta_VBOX_smooth(sim_truth(1):sim_truth(2)))
grid on
legend("mod","kin","wf","wf var","truth")



%Beta_VBOX_smooth(start_time_index:(start_time_index+length(out.Betay_mod.Data)-1),2)
[error_mean_mod,error_max_mod,time_at_max_mod,~] = errorCalc(out.Betay_mod.Data(sim(1):sim(2)),Beta_VBOX_smooth(sim_truth(1):sim_truth(2)));
disp(sprintf("MOD: MSE=%e; MAX=%e",error_mean_mod,error_max_mod))

[error_mean_kin,error_max_kin,time_at_max_kin,~] = errorCalc(out.Betay_kin.Data(sim(1):sim(2)),Beta_VBOX_smooth(sim_truth(1):sim_truth(2)));
disp(sprintf("KYN: MSE=%e; MAX=%e",error_mean_kin,error_max_kin))

[error_mean_wf,error_max_wf,time_at_max_wf,~] = errorCalc(out.Betay_wf.Data(sim(1):sim(2)),Beta_VBOX_smooth(sim_truth(1):sim_truth(2)));
disp(sprintf("WF: MSE=%e; MAX=%e",error_mean_wf,error_max_wf))

[error_mean_wf_var,error_max_wf_var,time_at_max_wf_var,~] = errorCalc(out.Betay_wf_var.Data(sim(1):sim(2)),Beta_VBOX_smooth(sim_truth(1):sim_truth(2)));
disp(sprintf("WF_VAR: MSE=%e; MAX=%e",error_mean_wf_var,error_max_wf_var))

%saving values to txt file
fileID = fopen('task1d.txt','a');
fprintf(fileID,"\n\n simulation: %s\n",sim_txt)
%fprintf(fileID,"time interval: [%f %f]\n",start_time_index)
fprintf(fileID,"Cf:%d    Cr:%d\n",Cf,Cr)
fprintf(fileID,"MOD: MSE=%e; MAX=%e\n",error_mean_mod,error_max_mod)
fprintf(fileID,"KYN: MSE=%e; MAX=%e\n",error_mean_kin,error_max_kin)
fprintf(fileID,"WF(T:%f): MSE=%e; MAX=%e\n",T,error_mean_wf,error_max_wf)
fprintf(fileID,"WF(T var): MSE=%e; MAX=%e\n",error_mean_wf_var,error_max_wf_var)
fclose(fileID);


