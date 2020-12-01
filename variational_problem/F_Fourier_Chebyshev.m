function [F] = F_Fourier_Chebyshev(c,ta,b,h,N,n,m)

%%% INPUTS
% c = (c_{ell,k})_{k=-N,...,N , ell = 0,...,n-1} : the Fourier-Chebyshev
%                    space-time coefficients of c (in k: Fourier and in ell: Chebyshev)
% ta = (ta_{ell,k})_{k=-N,...,N , ell = 0,...,n-1} : the Fourier-Chebyshev
%                    space-time coefficients of the approximate solution ta
% b = (b_k)_{k=-N,...,N} : row vector of Fourier coeff. of the initial condition
% angle : the angle parameter
% h : step size in time
% N : determines the # of Fourier coefficients which is 2*N+1
% n : # of Chebyshev coefficients

% theta = angle/180*pi;

omega = 2*pi;
ta = reshape(ta,n,2*N+1); ta = [ta;zeros(1,2*N+1)];

c = reshape(c,n,2*m+1); c = [c;zeros(1,2*m+1)];

%%% Beginning of the computation of the convolution

c_ext = [flipud(c(2:end,:));c];
ta_ext = [flipud(ta(2:end,:));ta];

ta_c = convtensor(ta_ext,[zeros(2*n+1,N-m),c_ext,zeros(2*n+1,N-m)]);

ta_c = ta_c(n+1:end,N+1-m:N+1+m);

%%% End of the computation of the convolution

F = zeros(n,2*m+1);

%%% Initial condition

for k = -m:m
    F(1,k+m+1) = c(1,k+m+1) + 2*sum( ((-1).^(1:n-1)').*c(2:n,k+m+1)) - b(k+m+1);
end

%%% The rest of the equations

for k = -m:m
    lambda_k = -(h/2)*(1i)*omega^2*k^2;
    ind = (1:n-1)';
    weight=ones(n-1,1); weight(n-1)=0;
    F(ind+1,k+m+1) = -lambda_k*c(ind,k+m+1)+2*ind.*c(ind+1,k+m+1)+lambda_k*c(ind+2,k+m+1)+h*(1i)*(weight.*ta_c(ind+2,k+m+1)-ta_c(ind,k+m+1));
end

F = reshape(F,(2*m+1)*n,1);

end