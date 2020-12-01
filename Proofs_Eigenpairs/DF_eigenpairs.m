function DF = DF_eigenpairs(x,v,theta)
%%%
%%% INPUT for F_eigenpairs: x = (lambda,a,b) in C^{2m+1}
%%% (1) lambda: the eigenvalue
%%% (2) a=(a0,a1,a2,...,a_{m-1}): vector of Fourier coefficients of the
%%% steady state
%%% (3) b=(b0,b1,b2,...,b_{m-1}): vector of Fourier coefficients of the
%%% eigenfunction
%%%

m = (length(x)-1)/2;

DF = zeros(2*m+1);

i0 = 1; i1 = (2:m+1); i2 = (m+2:2*m+1);

lambda = x(i0);
a = x(i1); %%% Fourier coefficients of the steady state
b = x(i2); %%% Fourier coefficients of the eigenfunction 
k = (0:m-1)';
omega = 2*pi;
mu = -k.^2*omega^2;

%%%%%%%%%%%%%
%%% eta_b %%%
%%%%%%%%%%%%%

DF(i0,i2) = v';

%%%%%%%%%%%%
%%% Df_a %%%
%%%%%%%%%%%%

a = [flip(a(2:end),1);a;zeros(m-1,1)];

Dfa_nonlinear = zeros(m);

Dfa_nonlinear(:,1) = 2*a(m:2*m-1);

for ell = 1:m-1
    Dfa_nonlinear(:,ell+1) = 2*(a(k-ell+m)+a(k+ell+m));
end 

Dfa = diag(mu) + Dfa_nonlinear;

DF(i1,i1) = Dfa;

%%%%%%%%%%%%%%%%%
%%% Dg_lambda %%%
%%%%%%%%%%%%%%%%%

DF(i2,i0) = -b;

%%%%%%%%%%%%
%%% Dg_a %%%
%%%%%%%%%%%%

b = [flip(b(2:end),1);b;zeros(m-1,1)];

Dga_nonlinear = zeros(m);

Dga_nonlinear(:,1) = 2*b(m:2*m-1);

for ell = 1:m-1
    Dga_nonlinear(:,ell+1) = 2*(b(k-ell+m)+b(k+ell+m));
end 

Dga = exp(1i*theta)*Dga_nonlinear;

DF(i2,i1) = Dga;

%%%%%%%%%%%%
%%% Dg_b %%%
%%%%%%%%%%%%

Dgb = exp(1i*theta)*(diag(mu) + Dfa_nonlinear) - lambda*eye(m);

DF(i2,i2) = Dgb;

end