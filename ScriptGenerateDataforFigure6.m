%This Matlab script generates the data required to plot SE for all
%different uplink receiver types in Figure 6.
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
nbrOfSetups = 500;

%Number of channel realizations per setup
nbrRealizations = 1000;

%H = sqrt(0.5)*(randn(N*L,K,nbrRealizations)+randn(N*L,K,nbrRealizations));
% Implemented for N = 1 case, perfect CSI

Q = diag(allocatedPowUEs);

%To avoid repeated reshapes of uplinkUEsPower.
powUEs = reshape(allocatedPowUEs,1,1,[]);
%R_tilde = zeros(size(R_tilde));

%
mseUEs_final_rls = zeros(nbrOfSetups,L);
mseUEs_final_sgd = zeros(nbrOfSetups,L);
mseUEs_final_ICC2020 = zeros(nbrOfSetups,L);
mseUEs_final_OSLP = zeros(nbrOfSetups,L);

for iter_nbrOfSetups = 1:nbrOfSetups
    
    disp(['Running Setup ' num2str(iter_nbrOfSetups) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate the simulation setup
    [gainOverNoisedB,R,pilotIndex,APpositions,UEpositions,distance_bw_UEs_APs] = generateSimulationSetup3(L,K,N,paramValues,channalDistribution,shadowFlag);
    
    %Compute channel estimates and Covariance matrices of estimates and its corresponding errors
    [H_hat,H,R_hat,R_tilde] = functionChannelEstimates3(R,nbrRealizations,L,K,N,tau_p,pilotIndex,allocatedPowUEs);
    
    mseUEs_sum_rls = zeros(K,L);
    mseUEs_sum_sgd = zeros(K,L);
    mseUEs_sum_ICC2020 = zeros(K,L);
    mseUEs_sum_OSLP = zeros(K,L);
    for iRealz = 1:nbrRealizations
        
        Sigma = sum(R_tilde(:,:,:,1).*powUEs,3) + eye(N);
        
        %RLS
        T_rls = eye(K);
        V_rls = 0;
        l = 1;
        [V_rls,T_rls] = functioncomputeV_RLS(H_hat(1:N,:,iRealz).',V_rls,T_rls,K,l);
        
        mseUEs_matrix_rls =  Q - (V_rls*H_hat(1:l*N,:,iRealz)*Q)' - (V_rls*H_hat(1:l*N,:,iRealz)*Q) + V_rls*( H_hat(1:l*N,:,iRealz)*diag(allocatedPowUEs)*(H_hat(1:l*N,:,iRealz)')+Sigma )*V_rls';
        mseUEs_sum_rls(:,l) = mseUEs_sum_rls(:,l) + diag(mseUEs_matrix_rls);
        
        %SGD
        V_sgd = 0;
        muSGD = 0.02;
        V_sgd = functioncomputeV_SGD(H_hat(1:N,:,iRealz).',V_sgd,muSGD,K,l);
        
        mseUEs_matrix_sgd =  Q - (V_sgd*H_hat(1:l*N,:,iRealz)*Q)' - (V_sgd*H_hat(1:l*N,:,iRealz)*Q) + V_sgd*( H_hat(1:l*N,:,iRealz)*diag(allocatedPowUEs)*(H_hat(1:l*N,:,iRealz)')+Sigma )*V_sgd';
        mseUEs_sum_sgd(:,l) = mseUEs_sum_sgd(:,l) + diag(mseUEs_matrix_sgd);
        
        %ICC2020
        V_ICC2020 = 0;
        vk_Hhat2_ICC2020 = 0;
        vk_Sigma2_vk_ICC2020 = 0;
        [V_ICC2020,vk_Hhat2_ICC2020,vk_Sigma2_vk_ICC2020] = functioncomputeV_ICC2020_2(H_hat(:,:,iRealz),R_tilde,V_ICC2020,vk_Hhat2_ICC2020,vk_Sigma2_vk_ICC2020,K,l,N,allocatedPowUEs);
        
        mseUEs_matrix_icc2020 = Q - (V_ICC2020*H_hat(1:l*N,:,iRealz)*Q)'- (V_ICC2020*H_hat(1:l*N,:,iRealz)*Q) + V_ICC2020*(H_hat(1:l*N,:,iRealz)*Q*(H_hat(1:l*N,:,iRealz)')+Sigma)*V_ICC2020';
        mseUEs_sum_ICC2020(:,l) = mseUEs_sum_ICC2020(:,l) + diag(mseUEs_matrix_icc2020);
        
        %OSLP
        P = Q;
        Hhat_1 = H_hat(1:N,:,iRealz);
        T_OSLP = (P*Hhat_1')/(Sigma + Hhat_1*P*Hhat_1');
        
        P = (eye(K) - T_OSLP*Hhat_1)*P;
        
        mseUEs_matrix_OSLP = P;
        mseUEs_sum_OSLP(:,l)= mseUEs_sum_OSLP(:,l) + diag(mseUEs_matrix_OSLP);
        
        SigmaL = Sigma; %K_L in paper
        for l = 2:L
            
            %Storing SigmaL Variable
            SigmaL = blkdiag(SigmaL,sum(R_tilde(:,:,:,l).*powUEs,3) + eye(N));
            
            %RLS
            [V_rls,T_rls] = functioncomputeV_RLS( H_hat((l-1)*N+1:l*N,:,iRealz).',V_rls,T_rls,K,l);
            
            mseUEs_matrix_rls =  Q - (V_rls*H_hat(1:l*N,:,iRealz)*Q)' - (V_rls*H_hat(1:l*N,:,iRealz)*Q) + V_rls*(H_hat(1:l*N,:,iRealz)*diag(allocatedPowUEs)*(H_hat(1:l*N,:,iRealz)')+SigmaL)*V_rls';
            mseUEs_sum_rls(:,l) = mseUEs_sum_rls(:,l) + diag(mseUEs_matrix_rls);
            
            %SGD
            V_sgd = functioncomputeV_SGD(H_hat((l-1)*N+1:l*N,:,iRealz).',V_sgd,muSGD,K,l);
            
            mseUEs_matrix_sgd =  Q - (V_sgd*H_hat(1:l*N,:,iRealz)*Q)' - (V_sgd*H_hat(1:l*N,:,iRealz)*Q) + V_sgd*( H_hat(1:l*N,:,iRealz)*diag(allocatedPowUEs)*(H_hat(1:l*N,:,iRealz)')+SigmaL )*V_sgd';
            mseUEs_sum_sgd(:,l) = mseUEs_sum_sgd(:,l) + diag(mseUEs_matrix_sgd);
            
            %ICC2020
            [V_ICC2020,vk_Hhat2_ICC2020,vk_Sigma2_vk_ICC2020] = functioncomputeV_ICC2020_2(H_hat(:,:,iRealz),R_tilde,V_ICC2020,vk_Hhat2_ICC2020,vk_Sigma2_vk_ICC2020,K,l,N,allocatedPowUEs);
            
            mseUEs_matrix_icc2020 = Q - (V_ICC2020*H_hat(1:l*N,:,iRealz)*Q)'- (V_ICC2020*H_hat(1:l*N,:,iRealz)*Q) + V_ICC2020*(H_hat(1:l*N,:,iRealz)*Q*(H_hat(1:l*N,:,iRealz)')+SigmaL)*V_ICC2020';
            mseUEs_sum_ICC2020(:,l) = mseUEs_sum_ICC2020(:,l) + diag(mseUEs_matrix_icc2020);
            
            %OSLP
            Hhat_2  = H_hat((l-1)*N+1:l*N,:,iRealz);
            Sigma_2 = sum(R_tilde(:,:,:,l).*powUEs,3); %Only for usage of OSLP
            
            T_OSLP = (P*Hhat_2')/(Sigma_2 + Hhat_2*P*Hhat_2');
            P = (eye(K) - T_OSLP*Hhat_2)*P;
            
            mseUEs_matrix_OSLP = P;
            mseUEs_sum_OSLP(:,l)= mseUEs_sum_OSLP(:,l) + diag(mseUEs_matrix_OSLP);
            
        end
        
    end
    
    %RLS
    mseUEs_avg_rls = mseUEs_sum_rls/nbrRealizations;
    mseUEs_avg_rls = real(mseUEs_avg_rls);
    
    mseRLS = sum(mseUEs_avg_rls,1)/trace(Q);
    
    %SGD
    mseUEs_avg_sgd = mseUEs_sum_sgd/nbrRealizations;
    mseUEs_avg_sgd = real(mseUEs_avg_sgd);
    
    mseSGD = sum(mseUEs_avg_sgd,1)/trace(Q);
    
    %ICC2020
    mseUEs_avg_ICC2020 = mseUEs_sum_ICC2020/nbrRealizations;
    mseUEs_avg_ICC2020 = real(mseUEs_avg_ICC2020);
    
    mseICC2020 = sum(mseUEs_avg_ICC2020,1)/trace(Q);
    
    
    %OSLP
    mseUEs_avg_OSLP = mseUEs_sum_OSLP/nbrRealizations;
    mseUEs_avg_OSLP = real(mseUEs_avg_OSLP);
    
    mseOSLP = sum(mseUEs_avg_OSLP,1)/trace(Q);
    
    % Storing values
    mseUEs_final_rls(iter_nbrOfSetups,:) = mseRLS;
    mseUEs_final_sgd(iter_nbrOfSetups,:) = mseSGD;
    mseUEs_final_ICC2020(iter_nbrOfSetups,:) = mseICC2020;
    mseUEs_final_OSLP(iter_nbrOfSetups,:) = mseOSLP;
    
end

%% Plots
mseUEs_final_rls = mean(mseUEs_final_rls,1);
mseUEs_final_sgd = mean(mseUEs_final_sgd,1);
mseUEs_final_ICC2020 = mean(mseUEs_final_ICC2020,1);
mseUEs_final_OSLP = mean(mseUEs_final_OSLP,1);

save('plotdataFigure6.mat');