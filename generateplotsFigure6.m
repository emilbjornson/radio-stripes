%This Matlab script generates the plots in Figure 6 of the work:
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

run('ScriptGenerateDataforFigure6.m');
%load('plotdataFigure6.mat');

MSE_OSLP_Final_50mW = mean(mseUEs_final_OSLP,1);
MSE_ICC2020_Final_50mW = mean(mseUEs_final_ICC2020,1);
MSE_RLS_Final_50mW = mean(mseUEs_final_rls,1);
MSE_SGD_Final_50mW = mean(mseUEs_final_sgd,1);

figure;
hold on; box on;
plot(1:L,10*log10(MSE_OSLP_Final_50mW),'-rp','linewidth',2,'MarkerSize',10);
plot(1:L,10*log10(MSE_ICC2020_Final_50mW),'-bs','linewidth',2,'MarkerSize',15);
plot(1:L,10*log10(MSE_RLS_Final_50mW),'-ko','linewidth',2,'MarkerSize',10);
legend('OSLP 50 mW','Algo. 2 50 mW','RLS 50 mW','FontSize',18,'location','best');
ylabel('Normalized MSE (dB)','FontSize',18);
xlabel('Number of APs L','FontSize',18);
%title(['Parameters: N = ',num2str(N),', L = ',num2str(L),', K = ',num2str(K),', and Radio Stripe Length = ',num2str(radioStripeLength)]);
grid on;