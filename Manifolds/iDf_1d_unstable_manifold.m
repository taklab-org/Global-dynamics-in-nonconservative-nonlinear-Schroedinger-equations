function Df = iDf_1d_unstable_manifold(p,par,r0)

%%% INPUTS
% p = (p_{k,m})_{k,m >= 0} : the Fourier-Taylor coefficients of P1
% par = (N,M,a,lambda,xi)
% N+1 is the Fourier projection
% M+1 is Taylor projection
% a : the fixed point (a vector of Fourier coefficients)
% lambda = the unstable eigenvalue
% xi : the unstable eigenvector

N = par(1); M = par(2);
theta = intval(par(2*N+6));

p = reshape(p,N+1,M+1);

tp = intval([p;zeros(N,M+1)]);

DN = intval(zeros((N+1)*(M+1)));

for m = 2:M
    m1 = (0:m);
    k = (0:N);
    for k1 = 0:N
        DN(k+1+m*(N+1),k1+1+m1*(N+1)) = 2*(1-isequal(k1,0)/2)*(tp(k+k1+1,m-m1+1)+tp(abs(k-k1)+1,m-m1+1));
    end
    
end

L = i_linear_operator_L(par,r0);

Df = L - exp(1i*theta)*DN;

end