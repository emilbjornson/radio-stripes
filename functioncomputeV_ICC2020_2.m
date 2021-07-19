function [Vbar,vk_Hhat2,vk_Sigma2_vk] = functioncomputeV_ICC2020_2(H_hat,R_tilde,Vbar,vk_Hhat2,vk_Sigma2_vk,K,l,N,allocatedPowUEs)
%This function computes N-LMMSE (Algorithm 2)receiver matrix
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


Q = diag(allocatedPowUEs);

powUEs = reshape(allocatedPowUEs,1,1,[]);

Sigma = sum(R_tilde(:,:,:,1).*powUEs,3) + eye(N);

if l==1
    
    Hhat1 = H_hat(1:N,:) ;
    
    V = (Q*Hhat1')/(Sigma + Hhat1*Q*Hhat1');% Is B as per general algorithm
    
    Vbar = V;
    
    vk_Hhat2     = V*Hhat1;
    vk_Sigma2_vk  = diag(V*Sigma*V');
    
else
    
    Hhat2 = [zeros(1,K);H_hat((l-1)*N+1:l*N,:)];
    
    Sigma2 = blkdiag(0,sum(R_tilde(:,:,:,l).*powUEs,3) + eye(N));
    
    V = zeros(K,N+1);
    for k = 1:K
        
        Hhat2(1,:)    = vk_Hhat2(k,:);
        Sigma2(1,1)   = vk_Sigma2_vk(k,1);
        
        % kth user combining vector
        vk2 = (Q(k,k)*Hhat2(:,k)')/(Sigma2+Hhat2*Q*Hhat2');
        
        V(k,:) = vk2;
        
        vk_Hhat2(k,:)     = vk2*Hhat2;
        vk_Sigma2_vk(k,1) = vk2*Sigma2*vk2';
        
    end
    
    A = diag(V(:,1));
    B = V(:,2:end);
    
    Vbar = [A*Vbar,B];
    
end


end