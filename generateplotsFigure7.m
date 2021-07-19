%This Matlab script generates the plots in Figure 7 of the work:
%
%Zakir Hussain Shaik, Emil Bjornson, and Erik G. Larsson,
%"MMSE-Optimal Sequential Processing for Cell-Free Massive MIMO With Radio
%Stripes," IEEE Transactions on Communications, To appear.
%
%Download article: https://arxiv.org/pdf/2012.13928.pdf
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.


clc;
clear;
close all;

maxL = 250;
L = 1:maxL;
N = 4;

K = 20;

tau_c = 200;
tau_p = 20;


L4_FH = 2*tau_c*N*L;

OSLP_FH = (2*K*(tau_c - tau_p) + K^2)*ones(1,maxL);
ICC2020_FH = (2*K*(tau_c - tau_p) + 2*K^2 + K)*ones(1,maxL);
SMR_FH = ( 2*K*(tau_c - tau_p) + K )*ones(1,maxL);

% Percentage of Fronthaul saved with L
OSLP_FH_saved = ( (L4_FH - OSLP_FH)./L4_FH )*100;
ICC2020_FH_saved = ( (L4_FH - ICC2020_FH)./L4_FH )*100;
SMR_FH_saved = ( (L4_FH - SMR_FH)./L4_FH )*100;

figure;
hold on;
plot(L,max(ICC2020_FH_saved,0),'-ro','linewidth',2);
plot(L,max(OSLP_FH_saved,0),'b','linewidth',2);
plot(L,max(SMR_FH_saved,0),'k--','linewidth',2);
legend('Algo. 2','OSLP/RLS','S-MR');
xlabel('Number of APs L','FontSize',18);
ylabel('Percentage of Fronthaul Saved','FontSize',18);
grid on;
box on;
