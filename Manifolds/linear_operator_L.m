function L = linear_operator_L(par)

N = par(1); M = par(2);
lambda = par(N+4);
theta = par(2*N+6);

L = zeros((N+1)*(M+1));
omega = 2*pi;
k = (0:N);
MU = diag(exp(1i*theta)*omega^2*k.^2);
Id = eye(N+1);

L(k+1,k+1) = Id; 
L(k+1+N+1,k+1+N+1) = Id;

for m=2:M
    L_m = lambda*m*Id + MU;
    L(k+1+m*(N+1),k+1+m*(N+1)) = L_m;
end

end

