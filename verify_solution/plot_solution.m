function plot_solution(a,tspan,index,N_pad,n_pad)
arguments
  a = []; tspan=[];index=1;
  N_pad = 100; n_pad = 100;
end
% if isintval(a)
%   a = mid(a);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input is the one-sided in Chebyshev %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chebfunpref.setDefaults('factory');
% index = 1;
[n,M] = size(a); %M=2*N+1
N = (M-1)/2;
N = N + N_pad; n = n + n_pad;
dx = 1/(2*N-1);
x = dx*(0:2*N);
a0 = [a(1,:);2*a(2:end,:);zeros(n_pad,M)];
a0 = [zeros(n,N_pad),a0,zeros(n,N_pad)];
if length(tspan)<2
  tspan=chebpts(n,[0,tspan]);
else
  tspan=chebpts(n,tspan);
end

if index==1
% Plot profile:
% figure
subplot(1,2,1);
surf(x,tspan,(2*N+1)*real(ifft(ifftshift(chebcoeffs2chebvals(a0),2).')).','EdgeColor','none')
hold on
xlabel('$x$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Re}(u)$','interpreter', 'latex')
title('Real part')
% figure
subplot(1,2,2);
surf(x,tspan,(2*N+1)*imag(ifft(ifftshift(chebcoeffs2chebvals(a0),2).')).','EdgeColor','none')
hold on
xlabel('$x$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Im}(u)$','interpreter', 'latex')
title('Imaginary part')
% 
elseif index==2
% figure
% Plot fourier modes:
k = (-N:N)';
subplot(1,2,1);
mesh(k,tspan,abs(real(chebcoeffs2chebvals(a0))))
% hold on
xlabel('$k$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Re}(u)$','interpreter', 'latex')
title('Real part')
set(gca, 'ZScale', 'log')
% figure
subplot(1,2,2);
mesh(k,tspan,abs(imag(chebcoeffs2chebvals(a0))))
% hold on
xlabel('$k$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Im}(u)$','interpreter', 'latex')
title('Imaginary part')
set(gca, 'ZScale', 'log')
elseif index==3
% norm plot
figure
LW = 'linewidth'; lw = 1.6;
y=(2*N+1)*ifft(ifftshift(chebcoeffs2chebvals(a0),2).');
plot(tspan,abs(max(y)),LW,lw);
% plot(tspan,max(abs(y)),LW,lw);
xlabel('$t$','interpreter','latex'), ylabel('$\|u(t)\|_{\infty}$','interpreter', 'latex')
title('The maximum norm of the solution.')
elseif index==4
% norm plot
figure
LW = 'linewidth'; lw = 1.6;
y=a0;
semilogy(abs(y),LW,lw);
xlabel('$\ell$','interpreter','latex'), ylabel('$|\tilde{a}_{\ell,k}|$','interpreter', 'latex')
title('Chebyshev coefficients.')
end