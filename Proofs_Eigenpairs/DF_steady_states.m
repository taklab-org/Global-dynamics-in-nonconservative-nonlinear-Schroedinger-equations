function DF = DF_steady_states(a,theta)

N = length(a)-1;
omega = 2*pi;

a = [flip(a(2:end),1);a;zeros(N,1)];

DF_nonlinear = zeros(N+1);

DF_nonlinear(:,1) = 2*a(N+1:2*N+1);

n = (0:N)';

for ell = 1:N
    DF_nonlinear(:,ell+1) = 2*(a(n-ell+N+1)+a(n+ell+N+1));
end 

DF = exp(1i*theta)*(diag(-n.^2*omega^2) + DF_nonlinear);

end

