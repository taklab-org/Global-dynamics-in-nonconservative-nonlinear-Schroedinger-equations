function Z0 = compute_Z0(B,m,n,nu)

% EEweights = nu(1).^(0:n-1)'*nu(2).^abs(-m:m);
% 
% [EEkt,EEkx] = ndgrid(0:n-1,-m:m);
% finite = (EEkt<n & abs(EEkx)<=m);
% weights = EEweights(finite);

nu = intval(nu);

weights = reshape(nu(1).^(0:n-1)'*nu(2).^abs(-m:m),n*(2*m+1),1);

Z0 = max((weights'*abs(B))./weights');

end

