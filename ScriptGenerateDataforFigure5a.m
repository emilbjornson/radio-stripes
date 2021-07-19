%This Matlab script generates the data required to plot SE for all
%different uplink receiver types in Figure 5a.
%
%This script was developed as a part of the paper:
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
clearvars;

%% System Setup
%Number of APs in the cell-free network
L = 24;

%Number of UEs in the network
K = 10;

%Number of antennas per AP
N = 1;

%% System Model Parameters-Specifications
%Communication bandwidth (in Hz)
B = 100e6;

%Length of the coherence block
tau_c = 200;

%Number of pilots per coherence block
tau_p = min(K,20);

%Uplink transmit power per UE (in W)
p=50/1000; %50 mW
allocatedPowUEs = p*ones(K,1);

%Noise figure (in dB)
noiseFigure = 9;

%Pathloss parameters
alpha = 3.67; % Pathloss exponent
PL_constantTerm = -30.5; % Median channel gain at a reference distance of 1Km in (dB)

%Standard deviation of the shadow fading
sigma_sf = 4;

%Decorrelation distance of the shadow fading
decorr = 9; % Distance beyond which UEs are decorrelated

%Height difference between an AP and a UE
distanceVertical = 5;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Angular standard deviation around the nominal angle (in degrees)
ASDdeg = 15;

%Size of the coverage area (as a square with wrap-around)
radioStripeLength = 500; %in meter, covering a square home 100 on each side

%Collecting all System Model Parameter Values
paramValues = {B,noiseFigure,tau_c,tau_p,allocatedPowUEs,alpha,PL_constantTerm,sigma_sf,decorr,distanceVertical,antennaSpacing,ASDdeg,radioStripeLength};

%% Monte Carlo Specifications

%Channel Distribution
%Enter accordingly as required: 'uncorrelated' , 'Gaussian' , 'Uniform' , 'Laplace'
channalDistribution = 'Gaussian';

%Modeling Shadowing then flag is 1 if not then it is 0
shadowFlag = 0;

%Number of Monte Carlo setups
nbrOfSetups = 100;

%Number of channel realizations per setup
numRealz = 200;

receiverType = 3; %receiverType = 1 for MR, receiverType = 2 for ZF and receiverType = 3 for N-LMMSE

SE_OSLP_Final = zeros(K,nbrOfSetups); %SINR of radio stripe
SE_RLS_Final = zeros(K,nbrOfSetups);
SE_SGD_Final = zeros(K,nbrOfSetups);
SE_ICC2020_Final = zeros(K,nbrOfSetups);
SE_SMR_Final = zeros(K,nbrOfSetups);
SE_CentMR_Final = zeros(K,nbrOfSetups);
SE_L4_Final = zeros(K,nbrOfSetups);

% Check of APs distance
if (radioStripeLength/L)<=50
    
    % msgbox('Warning ! APs are near than 50 m. Press any key to still continue');
    %pause;
    
end

muSGD = 0.02;
%% Monte Carlo over Setups
parfor iter_nbrOfSetups = 1:nbrOfSetups
    
    disp(['Running Setup ' num2str(iter_nbrOfSetups) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate the simulation setup
    [gainOverNoisedB,R,pilotIndex,APpositions,UEpositions,distance_bw_UEs_APs] = generateSimulationSetup3(L,K,N,paramValues,channalDistribution,shadowFlag);
    
    %Compute channel estimates and Covariance matrices of estimates and its corresponding errors
    [H_hat,H,R_hat,R_tilde] = functionChannelEstimates3(R,numRealz,L,K,N,tau_p,pilotIndex,allocatedPowUEs);
    
    %Calculate Spectral Efficiency of Radio Stripe in Uplink - using OSLP
    SE_OSLP = functioncomputeUplinkSE_OSLP(H_hat,R_tilde,tau_c,tau_p,numRealz,K,L,N,allocatedPowUEs);
    
    %Calculate Spectral Efficiency of Fully Centralized Level 4 Processing in Uplink
    SE_L4 = functioncomputeUplinkSE_L4(H_hat,R_tilde,tau_c,tau_p,numRealz,K,L,N,allocatedPowUEs);
    
    %Calculate Spectral Efficiency of Radio Stripe in Uplink - using S-MR
    SE_MR = functioncomputeUplinkSE_SMR(H_hat,R_tilde,tau_c,tau_p,numRealz,K,L,N,allocatedPowUEs);
    
    %Calculate Spectral Efficiency of Radio Stripe in Uplink - using S-MR
    SE_CentMR = functioncomputeUplinkSE_CentMR(H_hat,R_tilde,tau_c,tau_p,numRealz,K,L,N,allocatedPowUEs);
    
    %Calculate Spectral Efficiency of Radio Stripe in Uplink - using RLS
    SE_RLS = functioncomputeUplinkSE_RLS(H_hat,R_tilde,tau_c,tau_p,numRealz,K,L,N,allocatedPowUEs);
    
    %Calculate Spectral Efficiency of Radio Stripe in Uplink - using RLS
    SE_SGD = functioncomputeUplinkSE_SGD(H_hat,R_tilde,tau_c,tau_p,numRealz,K,L,N,allocatedPowUEs,muSGD);
    
    %Calculate Spectral Efficiency of Radio Stripe in Uplink - using ICC2020 Algorithm
    SE_ICC2020 = functioncomputeUplinkSE_ICC2020(H_hat,R_tilde,tau_c,tau_p,numRealz,K,L,N,allocatedPowUEs);
    
    SE_OSLP_Final(:,iter_nbrOfSetups) = SE_OSLP;
    SE_ICC2020_Final(:,iter_nbrOfSetups) = SE_ICC2020;
    SE_RLS_Final(:,iter_nbrOfSetups) = SE_RLS;
    SE_SGD_Final(:,iter_nbrOfSetups) = SE_SGD;
    SE_SMR_Final(:,iter_nbrOfSetups) = SE_MR;
    SE_CentMR_Final(:,iter_nbrOfSetups) = SE_CentMR;
    SE_L4_Final(:,iter_nbrOfSetups) = SE_L4;
    
end

save(['plotdataFigure5a_K',num2str(K),'.mat']);
