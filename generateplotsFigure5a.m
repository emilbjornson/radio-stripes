%This Matlab script generates the plots in Figure 5a of the work:
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
% To generate the results run the script: ScriptGenerateDataforFigure5a.m
% for K = 20, 24 which saves the data required to plot this


clc;
clear;
close all;


load('plotdataFigure5a_K20');
SE_OSLP_Final_K20P50mW = SE_OSLP_Final;
SE_ICC2020_Final_K20P50mW = SE_ICC2020_Final;
SE_RLS_Final_K20P50mW = SE_RLS_Final;
SE_SGD_Final_K20P50mW = SE_SGD_Final;
SE_SMR_Final_K20P50mW = SE_SMR_Final;
SE_ZF_Final_K20P50mW = SE_ZF_Final;

load('plotdataFigure5a_K24');
SE_OSLP_Final_K24P50mW = SE_OSLP_Final;
SE_ICC2020_Final_K24P50mW = SE_ICC2020_Final;
SE_RLS_Final_K24P50mW = SE_RLS_Final;
SE_SGD_Final_K24P50mW = SE_SGD_Final;
SE_SMR_Final_K24P50mW = SE_SMR_Final;
SE_ZF_Final_K24P50mW = SE_ZF_Final;

%%
close all;
figure;
hold on; grid on;
plot(sort(SE_RLS_Final_K24P50mW(:)),linspace(0,1,24*nbrOfSetups),'k','linewidth',2);
plot(sort(SE_ICC2020_Final_K24P50mW(:)),linspace(0,1,24*nbrOfSetups),'--b*','linewidth',2,'MarkerIndices',1:800:length(SE_ICC2020_Final_K24P50mW(:)),'MarkerSize',10);
plot(sort(SE_ICC2020_Final_K20P50mW(:)),linspace(0,1,20*nbrOfSetups),'-bd','linewidth',2,'MarkerIndices',1:800:length(SE_ICC2020_Final_K20P50mW(:)),'MarkerSize',10);

plot(sort(SE_OSLP_Final_K24P50mW(:)),linspace(0,1,24*nbrOfSetups),'--ro','linewidth',2,'MarkerIndices',1:800:length(SE_OSLP_Final_K24P50mW(:)),'MarkerSize',10);
plot(sort(SE_RLS_Final_K20P50mW(:)),linspace(0,1,20*nbrOfSetups),'--ks','linewidth',2,'MarkerIndices',1:800:length(SE_RLS_Final_K20P50mW(:)),'MarkerSize',10);
plot(sort(SE_OSLP_Final_K20P50mW(:)),linspace(0,1,20*nbrOfSetups),'--rp','linewidth',2,'MarkerIndices',1:800:length(SE_OSLP_Final_K20P50mW(:)),'MarkerSize',10);


ylabel('Cumulative Distribution','FontSize',18);
xlabel('Spectral efficiency (bit/s/Hz)','FontSize',18);
legend('RLS       K = 24','Algo. 2   K = 24','Algo. 2   K = 20','OSLP    K = 24','RLS       K = 20','OSLP    K = 20','location','best','FontSize',18);
box on;
