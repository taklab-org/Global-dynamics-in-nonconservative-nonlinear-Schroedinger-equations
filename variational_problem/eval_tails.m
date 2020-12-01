function [Fourier_tail_norm,Chebyshev_tail_norm] = eval_tails(c)

% Input c : n x (2*N+1) - Chebyshev (columns) x Fourier (rows) ceofficients 

c1 = norm(c(:,1:2),inf); c2 = norm(c(:,end-1:end),inf); % Fourier tail norms
c3 = norm(c(end-1:end,:),inf); % Chebyshev tail norm

Fourier_tail_norm = max([c1 c2]);
Chebyshev_tail_norm = c3;

end

