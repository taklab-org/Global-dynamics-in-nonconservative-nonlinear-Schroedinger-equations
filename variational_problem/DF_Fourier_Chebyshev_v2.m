function DF = DF_Fourier_Chebyshev_v2(ta,angle,h,N,n)

%%% INPUTS
% ta = (ta_{ell,k})_{k=-N,...,N , ell = 0,...,n-1} : the Fourier-Chebyshev
%                    space-time coefficients of the approximate solution ta
% angle : the angle parameter
% h : step size in time
% N : determines the # of Fourier coefficients which is 2*N+1
% n : # of Chebyshev coefficients

theta = angle/180*pi;

L = linear_operator_L(h,angle,N,n);

ta = reshape(ta,n,2*N+1); ta = [ta;zeros(n-1,2*N+1)];

DN = zeros((2*N+1)*n);

T = -diag(ones(1,n-1),-1) + diag(ones(1,n-1),1) ; T(1,2)=0;

for k = -N:N
    k1 = (-N:N);
    l = (0:n-1)';
    for l1 = 0:n-1
        S = (abs(k-k1)<=N);
        DN(l+1+(k+N)*n,l1+1+(k1(S)+N)*n) = (1-isequal(l1,0)/2)*h*exp(1i*theta)*T*(ta(abs(l-l1)+1,k-k1(S)+N+1)+ta(l+l1+1,k-k1(S)+N+1));
    end
end

DF = L + DN;

end
