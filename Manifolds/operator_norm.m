function norm_A = operator_norm(A,par,tail,nu)

N = par(1); M = par(2);
theta = intval(par(2*N+6));
omega = 2*intval('pi');
lambda_r = intval(real(par(N+4)));

k = (0:N);
m = (2:M)';

if tail > 0
    TAIL = [1;1;lambda_r*m+cos(theta)*(N+1)^2*omega^2].^(-1);
else
    TAIL = zeros(M+1,1);
end

alpha = zeros(M+1,1);

for m = 0:M
    alpha_m = 0;
    for j = 0:M
        A_jm = A(k+1+j*(N+1),k+1+m*(N+1));
        alpha_m = alpha_m + mnorma(A_jm,N,isequal(j,m)*TAIL(m+1),nu);
    end
    alpha(m+1) = sup(alpha_m);    
end

norm_A = max([max(alpha) sup(1/(lambda_r*(M+1)))]);

end