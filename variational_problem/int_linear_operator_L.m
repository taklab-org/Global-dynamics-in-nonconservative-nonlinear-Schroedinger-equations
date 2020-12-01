function [L] = int_linear_operator_L(h,N,n)

L = intval(zeros((2*N+1)*n));
ind = (1:n);
ipi = intval('pi');
omega = 2*ipi;
% theta = (angle/180)*ipi;
h = intval(h);

for k=-N:N
    lambda_k = -(h/2)*1i*omega^2*k^2;
    L_k = -lambda_k*diag(ones(1,n-1),-1) + diag(2*(ind-1)) + lambda_k*diag(ones(1,n-1),1) ;
    L_k(1)=1; L_k(1,2:2:end)=-2; L_k(1,3:2:end)=2;
    L(ind+(k+N)*n,ind+(k+N)*n) = L_k;
end

end

