function [tail_norm_fourier,tail_norm_chebyshev] = tail_norms(a,bar_k,bar_j)

a = reshape(a,bar_j,bar_k)';

tail_norm_fourier=norm(a(end-3:end,:));
tail_norm_chebyshev=norm(a(:,end-3:end));

display(['norm of Fourier tail = ',num2str(tail_norm_fourier),',  norm of Chebyshev tail = ',num2str(tail_norm_chebyshev)])

end

