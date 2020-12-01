function matrix_norm=mnorma(N,M,nu)

nu_power = [1 2*nu.^abs(1:N)];
matrix_norm=intval(max(mag(sum(abs(M).*repmat(nu_power,N+1,1)')./nu_power)));

end