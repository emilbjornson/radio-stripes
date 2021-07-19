function [H_hat,H,R_hat,R_tilde] = functionChannelEstimates3(R,numRealz,L,K,N,tau_p,pilotIndex,uplinkUEsPower)
%This function generate the channel realizations and estimates of these channels for all
%UEs in the entire network. The channels are modeled as correlated
%Rayleigh fading and the MMSE estimator is used.
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
%R                 = Matrix with dimension N x N x K x L where (:,:,k,l) is
%                    the spatial correlation matrix between AP l and UE k
%                    in setup n, normalized by the noise power
%numRealz          = Number of channel realizations
%L                 = Number of APs
%K                 = Number of UEs in the network
%N                 = Number of antennas per AP
%tau_p             = Number of orthogonal pilots
%pilotIndex        = Vector containing the pilot index assigned to each UE
%uplinkUEsPower    = Uplink transmit power per UE (same for everyone)
%
%OUTPUT:
%Hhat         = Matrix with dimension L*N x K x numRealz where
%               (:,k,n) is the estimated collective channel to UE k at
%               channel realization n.
%H            = Matrix with dimension L*N x K x numRealz with the
%               true channel realizations. The matrix is organized in the
%               same way as Hhat.
%R_hat        = Matrix with dimension N x N x K x L where (:,:,k,l) is the
%               spatial correlation matrix of the estimate between AP l and
%               UE k in setup n, normalized by the noise power
%R_tilde      = Matrix with dimension N x N x K x L where (:,:,k,l) is the
%               spatial correlation matrix of the estimatation error between AP l and
%               UE k in setup n, normalized by the noise power


%% Generate channel realizations

%Generate uncorrelated Rayleigh fading channel realizations
H = sqrt(0.5)*(randn(L*N,K,numRealz)+1i*randn(L*N,K,numRealz));

%Go through all channels and apply the spatial correlation matrices
for l = 1:L
    
    for k = 1:K
        
        %Apply correlation to the uncorrelated channel realizations
        Rsqrt = sqrtm(squeeze(R(:,:,k,l)));
        H((l-1)*N+1:l*N,k,:) = Rsqrt*squeeze(H((l-1)*N+1:l*N,k,:));
        
    end
    
end


%% Perform channel estimation

%Store identity matrix of size N x N
eyeN = eye(N);

%Generate realizations of normalized noise
n_tk = sqrt(0.5)*(randn(N,numRealz,L,tau_p) + 1i*randn(N,numRealz,L,tau_p));

%Prepare to store results
H_hat = zeros(L*N,K,numRealz);

R_hat = zeros(size(R));
R_tilde = zeros(size(R));


%Go through all APs
for l = 1:L
    
    %Go through all pilots
    for t = 1:tau_p
        
        %Compute processed pilot signal for all UEs that use pilot t
        yp = n_tk(:,:,l,t);
        Psi_tk = eyeN;
        for k = find(t==pilotIndex)'
            
            yp = yp + sqrt(tau_p*uplinkUEsPower(k))*reshape(H((l-1)*N+1:l*N,k,:),N,[]);
            Psi_tk = Psi_tk + tau_p*uplinkUEsPower(k)*R(:,:,k,l);
        end
        %Compute the matrix that is inverted in the MMSE estimator
        
        %Go through all UEs that use pilot t
        for k = find(t==pilotIndex)'
            
            R_PsiInv = R(:,:,k,l)/Psi_tk;
            
            %Compute the MMSE estimate
            H_hat((l-1)*N+1:l*N,k,:) = sqrt(tau_p)*sqrt(uplinkUEsPower(k))*R_PsiInv*yp;
            
            %Compute the spatial correlation matrix of the estimate
            R_hat(:,:,k,l) = uplinkUEsPower(k)*tau_p*R_PsiInv*R(:,:,k,l);
            R_tilde(:,:,k,l) = R(:,:,k,l) -  R_hat(:,:,k,l);
            
            
        end
        
    end
    
end

end
