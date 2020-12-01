function f = int_f_1d_unstable_manifold_ext(p,par,r0)

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
lambda = midrad(par(N+4),r0);
xi = par(N+5:2*N+5);
theta = intval(par(2*N+6));

omega = 2*intval('pi');

p = intval(reshape(p,N+1,M+1));

matrix_m = repmat((0:2*M),2*N+1,1);
matrix_k = repmat((0:2*N)',1,2*M+1);

p2_full = conv_TF(p,p);

p_ext = [[p intval(zeros(N+1,M))];intval(zeros(N,2*M+1))];

f = (lambda*matrix_m + exp(1i*theta)*matrix_k.^2*omega^2).*p_ext - exp(1i*theta)*p2_full;

f(:,1)=zeros(2*N+1,1); f(:,2)=zeros(2*N+1,1);

f(1:N+1,1) = p(:,1) - a;
f(1:N+1,2) = p(:,2) - xi;

f = reshape(f,(2*N+1)*(2*M+1),1);

end