%This Matlab script generates the plots in Figure 5b of the work:
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

run('ScriptGenerateDataforFigure5b');
%load('plotdataFigure5b.mat');

SE_OSLP_Final_1mW = SE_OSLP_Final;
SE_ICC2020_Final_1mW = SE_ICC2020_Final;
SE_RLS_Final_1mW = SE_RLS_Final;
SE_SGD_Final_1mW = SE_SGD_Final;
SE_SMR_Final_1mW = SE_SMR_Final;

%%
figure(3);
hold on; box on;
plot(sort(SE_SGD_Final_1mW(:)),linspace(0,1,K*nbrOfSetups),'--ks','linewidth',2,'MarkerIndices',1:400:length(SE_SGD_Final_1mW(:)),'MarkerSize',10);
plot(sort(SE_SMR_Final_1mW(:)),linspace(0,1,K*nbrOfSetups),'k','LineStyle','--','linewidth',2);
plot(sort(SE_RLS_Final_1mW(:)),linspace(0,1,K*nbrOfSetups),'-k*','linewidth',2,'MarkerIndices',1:400:length(SE_RLS_Final_1mW(:)),'MarkerSize',10);
plot(sort(SE_ICC2020_Final_1mW(:)),linspace(0,1,K*nbrOfSetups),'--bs','linewidth',2,'MarkerIndices',1:400:length(SE_ICC2020_Final(:)),'MarkerSize',10);
plot(sort(SE_OSLP_Final_1mW(:)),linspace(0,1,K*nbrOfSetups),'--ro','linewidth',2,'MarkerIndices',1:400:length(SE_OSLP_Final(:)),'MarkerSize',10);


legend('SGD     1 mW','S-MR    1 mW','RLS      1 mW','Algo. 2  1 mW','OSLP    1 mW','FontSize',18,'location','best');
ylabel('Cumulative Distribution','FontSize',18);
xlabel('Spectral efficiency (bit/s/Hz)','FontSize',18);
%title(['Parameters: N = ',num2str(N),', L = ',num2str(L),', K = ',num2str(K),', and Radio Stripe Length = ',num2str(radioStripeLength)]);
grid on;

%%
r_OSLP = sort(SE_OSLP_Final_1mW(:));
r_RLS = sort(SE_RLS_Final_1mW(:));
r_ICC2020 = sort(SE_ICC2020_Final_1mW(:));

percK10 = linspace(0,1,10*nbrOfSetups);

% 50% UEs
percentageUEs = 0.5;
[~,index] = min(abs(percK10 - percentageUEs));

gain_OSLPvRLS1 = r_OSLP(index) - r_RLS(index);
gain_ICC2020vRLS1 = r_ICC2020(index) - r_RLS(index);