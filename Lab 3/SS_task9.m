m = 22000; J = 700000;
c = 40000;
k = 600000;
k1 = k; k2= k;
L = 6;  L1= L; L2= L;
% cx=4e6;
% cz=24e3;
% cz=20.5e4;
% cx=50e5;
% cz = 250000;  %Robin and Gustaf's params
% cx = 700000;
cz=80e3;    %Lower force for impulse
cx=20e5;



% Wrong SS from before. This has been changed in the report
% A_ss_2dof_vh = [0  1 0 0;
%     -2*k/m, 0 , 0, 0;
%      0, 0, 0 ,1;
%      0, 0,-(2*k*L^2)/J,0];
% B_ss_2dof_vh = [0 0 0 0;
%                 k1/m k2/m  1/m  1/m;
%                  0     0     0     0;
%                  -(k1*L1)/J  k2*L2/J L1/J L2/J];
% C_ss_2dof_vh = [0 1 0 0;
%                 0 0 0 1];
% D_ss_2dof_vh = [0 0 0 0; 0 0 0 0];

A_ss_2dof_vh = [0 1 0 0;
    -2*k/m 0 0 0;
    0 0 0 1
    0 0 (-2*k*L^2)/J 0];

B_ss_2dof_vh = [0 0 0 0;
    k/m k/m -1/m -1/m;
    0 0 0 0;
    -k*L/J k*L/J L/J -L/J];

C_ss_2dof_vh = [0 1 0 0;
    0 0 0 1];

D_ss_2dof_vh = [0 0 0 0;
    0 0 0 0];
% out9 = sim("SS_task9",'StartTime','0','StopTime','3','FixedStep','0.1');
out9 = sim("SS_task9_2021",'StartTime','0','StopTime','20','FixedStep','0.01');