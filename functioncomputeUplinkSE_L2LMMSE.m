function SE = functioncomputeUplinkSE_L2LMMSE(H_hat,R_tilde,tau_c,tau_p,numRealz,K,L,N,allocatedPowUEs)
%This function computes achievable SE of a radio stripes network with
%Local LMMSE receiver
%
%This function was developed as a part of the paper:
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
%
%INPUT:
%Hhat                   = Matrix with dimension L*N x K x numRealz where
%                         (:,k,n) is the estimated collective channel to UE k at
%                          channel realization n.
%R_tilde                = Matrix with dimension N x N x K x L where (:,:,k,l) is the
%                         spatial correlation matrix of the estimatation error between AP l and
%                         UE k in setup n, normalized by the noise power
%tau_c                  = Length of the coherence block
%tau_p                  = Number of channel uses for piloting
%numRealz               = Number of channel realizations
%K                      = Number of UEs in the network
%L                      = Number of APs for the Radio-Stripe Network
%N                      = Number of antennas per AP
%allocated powers       = Power allocated to UEs, K x 1 vector
%
%OUTPUT:
%SE                     = Spectral Efficiency vector of K UEs.
%
%Generates SE data (K x 1) vector where SE(k) is achievable SE of UE k
%with centralized Local LMMSE prcoessing


% Variable to store the final result
SE = zeros(K,1);

% Transmit vector covariance matrix
Q = diag(allocatedPowUEs);

% Reshaping power coefficient vector to 3rd dimension
powUEs = reshape(allocatedPowUEs,1,1,[]);
% To avoid repeated reshapes of allocatedPowUEs.
% This reshape to 3rd dimension is done because: R_tilde's variable 3rd
% dimension is UEs and SE expression requires weighted sum of
% error covariance matrices over all UEs where weights are power
% coefficients. Hence, it would be easy to obtain the summation as below:

% Variable to store summation of covariance matrices (over UEs) of all APs as block
% diagonal. Because here we utilize centralized network SE expression.
Sigma = sum(R_tilde(:,:,:,1).*powUEs,3) + eye(N);

% Iterate over channel realizations
for iRealz = 1:numRealz
    
    Hhat_1 = H_hat(1:N,:,iRealz); % Collecting channel estimate of AP1
    
    B = (1/L)*(Q*Hhat_1')/(Sigma + Hhat_1*Q*Hhat_1'); % B1 , eq. (20)
    
    V = B; % For consistency in all code files, V is used instead of B
    
    Sigma2 = Sigma; %K_L in paper
    for l = 2:L
        
        Hhat_2  = H_hat((l-1)*N+1:l*N,:,iRealz);
        Sigma = sum(R_tilde(:,:,:,l).*powUEs,3) + eye(N);
        
        A = eye(K);
        B = (1/L)*(Q*Hhat_2')/(Sigma + Hhat_2*Q*Hhat_2');
        
        V = [A*V,B]; % Eq. (20)
        
        %Storing Sigma2 Variable
        Sigma2 = blkdiag(Sigma2,Sigma); % K_L
        
    end
    
    for k = 1:K
        
        hk_hat = H_hat(:,k,iRealz); % UE k estimate channel vector in eq. (23)
        vk     = V(k,:)'; % UE k combining vector in eq. (23) (note: bk = vk)
        
        % Computing Numerator and Denominator of eq. (23)
        sinr_numer = ( allocatedPowUEs(k)*abs(vk'*hk_hat)^2 );
        sinr_denom = vk'*(H_hat(:,:,iRealz)*diag(allocatedPowUEs)*H_hat(:,:,iRealz)')*vk - sinr_numer + vk'*Sigma2*vk;
        
        % Sum rate over all realizations
        SE(k,1) = SE(k,1) + log2(1 + real(sinr_numer/sinr_denom) );
        
    end
    
end

SE = (1 - tau_p/tau_c)*SE/numRealz; % Average rate

end