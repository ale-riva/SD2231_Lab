%% This is a Matlab file for designing H_infinity controller (assignment 3
%% of SD2231)
s=tf('s');

% systme parameters
m=1.15*22000;   %kg
j=1.15*700e3;   %kgm^2
c=0.85*40e3;    %Ns/m
k=0.85*2*300e3; %N/m
L=6;       %m

%% State space model for skyhook contorl
Ask=[0 1 0 0
    -2*k/m 0 0 0
    0 0 0 1
    0 0 -2*k*L^2/j 0];
Bsk=[0 0 0 0
    k/m k/m -1/m -1/m
    0 0 0 0
    -L*k/j L*k/j L/j -L/j];
Csk=[0 1 0 0
    0 0 0 1];
Dsk=zeros(2,4);

%% H_inf using linmod syntax

%state space: The same as skyhook

      
%Weighting functions

%For penalizing actuator force
Wa1=(0.00175*s+1)/(0.00025*s+1);
Wa2=Wa1;

%For penalizing bounce and pitch motions
eps=1;
wn = eig(Ask);         %Outputs the eigenvalues/natural frequencies
wnb=7.3855;            %Find the right equation or value for wnb
wnchi=7.8558;          %Find the right equation or value for wnchi
s1b=-eps+1i*sqrt(wnb^2-eps^2);
s2b=-eps-1i*sqrt(wnb^2-eps^2);
s1chi=-eps+1i*sqrt(wnchi^2-eps^2);
s2chi=-eps-1i*sqrt(wnchi^2-eps^2);
% kb=input('Enter the gain for Wb = '); 
% kchi=input('Enter the gain for Wchi = ');
kb = 7000;
kchi = 5000;
Wb=(kb*s1b*s2b)/((s-s1b)*(s-s2b));
Wchi=(kchi*s1chi*s2chi)/((s-s1chi)*(s-s2chi));

%Extracting the extended model
[A_Pe,B_Pe,C_Pe,D_Pe] = linmod('Extended_model_s');% state space parameters of the extended system: Pe
Pe=ss(A_Pe,B_Pe,C_Pe,D_Pe);

%Calculating the controller
ncont = 2;%Number of control inputs
nmeas = 2;%Number of measured outputs provided to the controller
Pe=minreal(Pe);%This syntax cancels pole-zero pairs in transfer
%functions. The output system has minimal order and the same response
%characteristics as the original model.
%[K,Pec,gamma,info]=hinfsyn(Pe,nmeas,ncont,'method','lmi'); % for working with the error
%[Ainf, Binf, Cinf, Dinf]=ssdata(K);
%%
out11_3_h = sim("H_inf_2021",'StartTime','0','StopTime','20','FixedStep','0.01');

%%
%Now use the controller K in your simulation
figure
bodemag(Wa1, Wa2, Wb, Wchi)
xline(7.3855)
hold on
xline(7.8558)
legend("Wa1","Wa2","Wb","Wchi")
