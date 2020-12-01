function nu_norm = compute_nu_norm(c,nu)

[n,tN] = size(c); N = (tN-1)/2;

nu1 = nu(1); nu2 = nu(2);

%ext_c = [flipud(c(2:n,:));c];

ell = abs(0:n-1)'; k = abs(-N:N); 

weights = (nu1.^ell)*(nu2.^k);

nu_norm = sum(sum(abs(c).*weights));

% [n,tN] = size(c); N = (tN-1)/2;
% 
% nu1 = nu(1); nu2 = nu(2);
% 
% ext_c = [flipud(c(2:n,:));c];
% 
% ell = abs(-n+1:n-1)'; k = abs(-N:N); 
% 
% weights = (nu1.^ell)*(nu2.^k);
% 
% nu_norm = sum(sum(abs(ext_c).*weights));

end

