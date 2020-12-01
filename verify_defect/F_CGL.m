function [F] = F_CGL(a,angle,tmax,N,n)

a = reshape(a,n,2*N+1);

theta = angle/180*pi;

[~,n2] = size(a);
N = (n2-1)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Computation of the coefficients of d/dt a_k(t) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

da = (2/tmax)*ChebDerCoeffs_fourier(a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beginning of the computation of the convolution %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = [flipud(a(2:end,:));a(1,:);a(2:end,:)];
a2 = convtensor(a,a);
a = a(n:end,:);
a2 = a2(n:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% End of the computation of the convolution %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = -N:N;
F = da + exp(1i*theta)*((4*pi^2)*(k.^2).*a-a2);
F = reshape(F,(2*N+1)*n,1);

end

