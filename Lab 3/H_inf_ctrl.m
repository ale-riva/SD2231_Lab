%% This is a Matlab file for designing H_infinity controller (assignment 3
%% of SD2231)
clear all
close all
clc

s=tf('s');

% systme parameters
m=22000;   %kg
j=700e3;   %kgm^2
c=40e3;    %Ns/m
k=2*300e3; %N/m
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

%bode plot of Weight functions of forces
w=logspace(0,6,2*100040);
[A_wa1,phi_wa1] = bode(Wa1,w);
figure
subplot(2,1,1)
loglog(w,A_wa1(:))
grid on
subplot(2,1,2)
semilogx(w,phi_wa1(:))
grid on


%%
%frequency we should penalize the most-- the one corresponding to the peaks
w_nat = abs(eig(Ask));



%For penalizing bounce and pitch motions
eps=1;
wnb=w_nat(1);            %Find the right equation or value for wnb
wnchi=w_nat(3);          %Find the right equation or value for wnchi
s1b=-eps+1i*sqrt(wnb^2-eps^2);
s2b=-eps-1i*sqrt(wnb^2-eps^2);
s1chi=-eps+1i*sqrt(wnchi^2-eps^2);
s2chi=-eps-1i*sqrt(wnchi^2-eps^2);
%kb=input('Enter the gain for Wb = '); 
%kchi=input('Enter the gain for Wchi = ');

kb=7e3;
kchi = 20e3;
Wb=(kb*s1b*s2b)/((s-s1b)*(s-s2b));
Wchi=(kchi*s1chi*s2chi)/((s-s1chi)*(s-s2chi));

% w=logspace(-2,4,2*100040);
% [A_wb,phi_wb] = bode(Wb,w);
% figure
% subplot(2,1,1)
% loglog(w,A_wb(:))
% grid on
% subplot(2,1,2)
% semilogx(w,phi_wb(:))
% grid on
% sgtitle("Wb")
% 
% [A_wchi,phi_wchi] = bode(Wchi,w);
% figure
% subplot(2,1,1)
% loglog(w,A_wchi(:))
% grid on
% subplot(2,1,2)
% semilogx(w,phi_wchi(:))
% grid on
% sgtitle("Wchi")
%%
%Extracting the extended model
[A_Pe,B_Pe,C_Pe,D_Pe] = linmod('Extended_model');% state space parameters of the extended system: Pe
Pe=ss(A_Pe,B_Pe,C_Pe,D_Pe);
%%
%Calculating the controller
ncont = 2;%Number of control inputs
nmeas = 2;%Number of measured outputs provided to the controller
Pe=minreal(Pe);%This syntax cancels pole-zero pairs in transfer
%functions. The output system has minimal order and the same response
%characteristics as the original model.
[K,Pec,gamma,info]=hinfsyn(Pe,nmeas,ncont,'method','lmi'); % for working with the error
[Ainf, Binf, Cinf, Dinf]=ssdata(K);

%Now use the controller K in your simulation