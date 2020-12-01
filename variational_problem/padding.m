function [a,b,bar_k,bar_j]=padding(a,b,bar_k,bar_j,d1,d2)

%%% bar_k : size of the Fourier expansion
%%% bar_j : size of the Chebyshev expansion

%%% d1: extra_dimensions in Fourier (can be negative)
%%% d2: extra_dimensions in Chebyshev (can be negative)
%%% (we assume that d1 and d2 have the same sign)

a = reshape(a,bar_j,bar_k)';

a_ext = zeros(bar_k+d1,bar_j+d2);
b_ext = zeros(bar_k+d1,1);

if d1>=0,
    a_ext(1:bar_k,1:bar_j) = a(1:bar_k,1:bar_j);
    b_ext(1:bar_k,1) = b(1:bar_k,1);
else
    a_ext=a(1:bar_k+d1,1:bar_j+d2);    
    b_ext=b(1:bar_k+d1,1);
end

a = a_ext;
b = b_ext;
bar_k=bar_k+d1;
bar_j=bar_j+d2;

a=reshape(a',bar_k*bar_j,1);

end