function DF = finite_diff_F_Fourier_Chebyshev_backward(c,ta,b,angle,h,N,n,m)

delta = 1e-8;
m1 = length(c);
E = eye(m1);
DF = zeros(m1);

for j=1:m1
    ch = c + delta*E(:,j);
    DF(:,j)=(F_Fourier_Chebyshev_backward(ch,ta,b,angle,h,N,n,m)-F_Fourier_Chebyshev_backward(c,ta,b,angle,h,N,n,m))/delta;
end
end