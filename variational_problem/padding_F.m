function [F,bar_k,bar_j]=padding_F(F,bar_k,bar_j,d1,d2)

%%% bar_k : size of the Fourier expansion
%%% bar_j : size of the Chebyshev expansion

%%% d1: extra_dimensions in Fourier (can be negative)
%%% d2: extra_dimensions in Chebyshev (can be negative)
%%% (we assume that d1 and d2 have the same sign)

F = reshape(F,bar_j,bar_k)';

F_ext = zeros(bar_k+d1,bar_j+d2);

if d1>=0,
    F_ext(1:bar_k,1:bar_j) = F(1:bar_k,1:bar_j);
else
    F_ext=F(1:bar_k+d1,1:bar_j+d2);    
end

F = F_ext;
bar_k=bar_k+d1;
bar_j=bar_j+d2;

F=reshape(F',bar_k*bar_j,1);

end