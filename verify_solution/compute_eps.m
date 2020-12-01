function eps=compute_eps(a,v_hat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input is the one-sided in Chebyshev %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chebfunpref.setDefaults('factory');
% index = 1;
[n,M] = size(a); %M=2*N+1
N = (M-1)/2;
% dx = 1/(2*N-1);
% x = dx*(0:2*N);
a0 = [a(1,:);2*a(2:end,:)];
%
ell = 0:n-1;
v = (-1).^ell;
% v_hat = zeros(1,2*N+1);
% v_hat(N+1)=50; v_hat(N)=-25; v_hat(N+2)=-25;
% 
eps = norm(v*a0-v_hat,1);