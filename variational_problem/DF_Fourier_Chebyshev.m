function DF = DF_Fourier_Chebyshev(ta,h,N,n,core)

%%% INPUTS
% ta = (ta_{ell,k})_{k=-N,...,N , ell = 0,...,n-1} : the Fourier-Chebyshev
%                    space-time coefficients of the approximate solution ta
% angle : the angle parameter
% h : step size in time
% N : determines the # of Fourier coefficients which is 2*N+1
% n : # of Chebyshev coefficients

% theta = angle/180*pi;

ta = reshape(ta,n,2*N+1); ta = [ta;zeros(n-1,2*N+1)];

% N = core;

L = linear_operator_L(h,core,n);

DN = zeros((2*core+1)*n);

T = -diag(ones(1,n-1),-1) + diag(ones(1,n-1),1) ; T(1,2)=0;

for k = -core:core
    for k1 = -core:core
        l = (0:n-1)';
        l1 = 0;
        if abs(k-k1)<=N
            DN(l+1+(k+core)*n,l1+1+(k1+core)*n) = h*(1i)*T*ta(abs(l-l1)+1,k-k1+N+1);
        end
        for l1 = 1:n-1         
            if abs(k-k1)<=N
                DN(l+1+(k+core)*n,l1+1+(k1+core)*n) = h*(1i)*T*(ta(abs(l-l1)+1,k-k1+N+1)+ta(l+l1+1,k-k1+N+1));
            end
        end
    end
end

DF = L + DN;

end
