function matrix_norm = mnorma(B,N,tail,nu)

nu_power = [1 2*nu.^abs(1:N)];
C1 = max(mag(sum(abs(B).*repmat(nu_power,N+1,1)')./nu_power));
C2 = tail;
matrix_norm = intval(max([C1 C2]));

end