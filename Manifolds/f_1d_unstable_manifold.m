function f = f_1d_unstable_manifold(p,par)

%%% INPUTS
% p = (p_{k,m})_{k,m >= 0} : the Fourier-Taylor coefficients of P1
% par = (N,M,a,lambda,xi)
% N+1 is the Fourier projection
% M+1 is Taylor projection
% a : the fixed point (a vector of Fourier coefficients)
% lambda = the unstable eigenvalue
% xi : the unstable eigenvector

N = par(1); M = par(2);
a = par(3:N+3); 
lambda = par(N+4);
xi = par(N+5:2*N+5);
theta = par(2*N+6);

omega = 2*pi;

p = reshape(p,N+1,M+1);

matrix_m = repmat((0:M),N+1,1);
matrix_k = repmat((0:N)',1,M+1);

[~,p2] = conv_TF(p,p);

f = (lambda*matrix_m + exp(1i*theta)*matrix_k.^2*omega^2).*p - exp(1i*theta)*p2;

f(:,1) = p(:,1) - a;
f(:,2) = p(:,2) - xi;

f = reshape(f,(N+1)*(M+1),1); 

end