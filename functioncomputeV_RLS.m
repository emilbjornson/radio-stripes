function [V,T] = functioncomputeV_RLS(h,V,T,K,l)
%This function computes RLS receiver matrix
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


z = T*conj(h);
a = 1/(1 + (h.')*z);

T = T - a*(z*z');

A = (eye(K) - a*z*transpose(h));
B = a*z;

if l==1
    
    V = B;
    
else
    V = [A*V,B];
end

end