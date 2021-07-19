%This Matlab script generates the plots in Figure 4a of the work:
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

run('ScriptGenerateDataforFigure4a');
%load('plotdataFigure4a');

figure;
hold on; box on;grid on;
plot(sort(SE_OSLP_Final(:)),linspace(0,1,K*nbrOfSetups),'k*','linewidth',2,'MarkerIndices',1:300:length(SE_OSLP_Final(:)),'MarkerSize',15);
plot(sort(SE_L2LLMSE_Final(:)),linspace(0,1,K*nbrOfSetups),'-x','linewidth',2,'MarkerIndices',1:300:length(SE_L2LLMSE_Final(:)),'MarkerSize',15);
plot(sort(SE_ICC2020_Final(:)),linspace(0,1,K*nbrOfSetups),'-ks','linewidth',2,'MarkerIndices',1:300:length(SE_ICC2020_Final(:)),'MarkerSize',10);
plot(sort(SE_SMR_Final(:)),linspace(0,1,K*nbrOfSetups),'ko','linewidth',2,'MarkerIndices',1:300:length(SE_SMR_Final(:)),'MarkerSize',10);
plot(sort(SE_CentMR_Final(:)),linspace(0,1,K*nbrOfSetups),'--r','linewidth',2);
plot(sort(SE_L4_Final(:)),linspace(0,1,K*nbrOfSetups),'b','linewidth',2);
legend('OSLP','Local LMMSE','Algo. 2','S-MR','Cent MR','Cent LMMSE','FontSize',18,'location','best');
ylabel('Cumulative Distribution','FontSize',18);
xlabel('Spectral efficiency (bit/s/Hz)','FontSize',18);