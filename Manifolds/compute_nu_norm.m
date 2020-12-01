function nu_norm = compute_nu_norm(p,par,nu)

N = par(1); M = par(2);

p = reshape(p,N+1,M+1);

k = abs(0:N)'; 

weights = (nu.^k)*ones(1,M+1);

nu_norm = sum(sum(abs(p).*weights));

end

