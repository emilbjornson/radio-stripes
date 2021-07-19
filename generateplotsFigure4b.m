%This Matlab script generates the plots in Figure 4b of the work:
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

%%%%% Important Note
% To generate the results run the script: ScriptGenerateDataforFigure4b.m
% for K = 5, 10, 20, 30, 40, 48 which saves the data required to plot this


clc;
clear;
close all;


L = 12;
KK = [5,10,20,30,40,48];

avgR_MR = zeros(length(KK),1);
avgR_ICC2020 = zeros(length(KK),1);
avgR_L2LLMSE = zeros(length(KK),1);
avgR_L4 = zeros(length(KK),1);
avgR_OSLP = zeros(length(KK),1);
avgR_SMR = zeros(length(KK),1);
avgR_ZF = zeros(length(KK),1);

for k = 1:length(KK)
    
    load(['plotdataFigure4b_K',num2str(KK),'.mat']);
    
    avgR_MR(k,1) = mean(SE_CentMR_Final(:));
    avgR_ICC2020(k,1) = mean(SE_ICC2020_Final(:));
    avgR_L2LLMSE(k,1) = mean(SE_L2LLMSE_Final(:));
    avgR_L4(k,1) = mean(SE_L4_Final(:));
    avgR_OSLP(k,1) = mean(SE_OSLP_Final(:));
    avgR_SMR(k,1) = mean(SE_SMR_Final(:));
    avgR_ZF(k,1) = mean(SE_ZF_Final(:));
    
    
end

%
plot(KK,avgR_L4,'-o','linewidth',2); hold on;
plot(KK,avgR_OSLP,'kd','MarkerSize',15);
plot(KK,avgR_ZF,'--+','linewidth',2,'MarkerSize',10);
plot(KK,avgR_ICC2020,'-x','linewidth',2,'MarkerSize',14);
plot(KK,avgR_L2LLMSE,'-s','linewidth',2,'MarkerSize',14);
plot(KK,avgR_SMR,'-*','linewidth',2,'MarkerSize',14);
legend('Cent LMMSE','OSLP','Cent ZF','Algo. 2','Local LMMSE','S-MR','FontSize',18,'location','best');
xlabel('Number of UEs (K)');
ylabel('Average rate (bit/s/Hz)');
grid on;
xlim([min(KK),max(KK)]);