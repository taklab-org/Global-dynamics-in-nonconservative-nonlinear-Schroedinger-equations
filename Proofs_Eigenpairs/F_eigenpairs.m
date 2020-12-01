function F = F_eigenpairs(x,v,theta)
%%%
%%% INPUT for F_eigenpairs: x = (lambda,a,b) in C^{2m+1}
%%% (1) lambda: the eigenvalue
%%% (2) a=(a0,a1,a2,...,a_{m-1}): vector of Fourier coefficients of the
%%% steady state
%%% (3) b=(b0,b1,b2,...,b_{m-1}): vector of Fourier coefficients of the
%%% eigenfunction
%%%
%%% OUTPUT for F_eigenpair: F=(eta,f,g) in C^{2m+1}
%%% (1) eta : the phase condition for the eigenvector b
%%% (2) f=(f0,f1,f2,...,f_{m-1})
%%% (3) g=(g0,g1,g2,...,g_{m-1})

m = (length(x)-1)/2;

lambda = x(1);
a = x(2:m+1); %%% Fourier coefficients of the steady state
b = x(m+2:2*m+1); %%% Fourier coefficients of the eigenfunction 
k = (0:m-1)';
omega = 2*pi;
mu = -k.^2*omega^2;

a2 = quadraticFFT(a,a);
ab = quadraticFFT(a,b);

eta = v'*b - 1;
f = mu.*a + a2 ;
g = exp(1i*theta)*(mu.*b + 2*ab) - lambda*b ;

F = [eta;f;g];

end