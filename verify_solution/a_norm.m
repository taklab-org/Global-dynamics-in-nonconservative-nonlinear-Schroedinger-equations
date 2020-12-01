function norma = a_norm(a)

% Input: one-sided Chebyshev coefficients, (a_{ell,k})_{ell=0,...,n-1, |k|<=N}
% Output: X-norm of a, i.e., \sup_t \|a(t)\|_{\ell^1}
a = [a(1,:);2*a(2:end,:)]; % Back to Chebyshev coefficient

norma = sum(sum(abs(a)));
