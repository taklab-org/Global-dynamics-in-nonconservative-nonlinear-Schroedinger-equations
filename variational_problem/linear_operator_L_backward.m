function [L] = linear_operator_L_backward(h,N,n)

L = zeros((2*N+1)*n);
ind = (1:n);
omega = 2*pi;
% theta = angle/180*pi;

for k=-N:N
    lambda_k = -(h/2)*(1i)*omega^2*k^2;
    L_k = diag(2*(ind-1)) - (-lambda_k*diag(ones(1,n-1),-1) + lambda_k*diag(ones(1,n-1),1));
    L_k(1)=1; L_k(1,2:2:end)=-2; L_k(1,3:2:end)=2;
    L(ind+(k+N)*n,ind+(k+N)*n) = L_k;
end

end

