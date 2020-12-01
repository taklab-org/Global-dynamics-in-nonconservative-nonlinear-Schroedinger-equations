function [tspan,y2]=plot_solution2(a,tspan,index,phi0)
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
dx = 1/(2*N-1);
x = dx*(0:2*N);
a0 = [a(1,:);2*a(2:end,:)];
if length(tspan)<2
  tspan=chebpts(n,[0,tspan]);
else
  tspan=chebpts(n,tspan);
end
if index==1
% Plot profile:
% figure
subplot(1,2,1);
mesh(x,tspan,(2*N+1)*real(ifft(ifftshift(chebcoeffs2chebvals(a0),2).')).')
hold on
xlabel('$x$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Re}(u)$','interpreter', 'latex')
title('Real part')
% figure
subplot(1,2,2);
mesh(x,tspan,(2*N+1)*imag(ifft(ifftshift(chebcoeffs2chebvals(a0),2).')).')
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
hold on
xlabel('$k$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Re}(u)$','interpreter', 'latex')
title('Real part')
set(gca, 'ZScale', 'log')
% figure
subplot(1,2,2);
mesh(k,tspan,abs(imag(chebcoeffs2chebvals(a0))))
hold on
xlabel('$k$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Im}(u)$','interpreter', 'latex')
title('Imaginary part')
set(gca, 'ZScale', 'log')
elseif index==3
% norm plot
figure
LW = 'linewidth'; lw = 1.6;
y=ifft(ifftshift(chebcoeffs2chebvals(a0),2).');
plot(tspan,abs(max(y)),LW,lw);
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
elseif index==5
ta=chebcoeffs2chebvals(a0);
y1=abs(ta(:,N+1));
y2=abs(ta(:,N+2));
y3=abs(ta(:,N));
y4=sum(abs(ta(:,1:N)),2)+sum(abs(ta(:,N+2:end)),2);
semilogy(tspan,abs([y4./y1.^2]));
xlabel('$t$','interpreter','latex'), ylabel('$\|u(t)\|_{\infty}$','interpreter', 'latex')
elseif index==6
ta=chebcoeffs2chebvals(a0);
y1=ta(:,N+1);% a0(t)
% phi=y1(1);
% ep=@(t) -phi0./(-1+1i*phi0*t);
% y4=sum(abs(ta(:,1:N)),2)+sum(abs(ta(:,N+2:end)),2);
% y2=abs(y1./ep(tspan)-1);

% loglog(tspan,abs(ta(:,N-4:N+1)));
plot(tspan,abs(ta(:,N+1)),'b-');
plot(tspan,abs(ta(:,N  )),'r-');
plot(tspan,abs(ta(:,N-1)),'g-');
% loglog(tspan,abs(ta(:,N-2)),'k-');
% loglog(tspan,abs(ta(:,N-3)),'m-');
% loglog(tspan,abs(ta(:,N-4)),'c-');
% plot(tspan,abs(ep(tspan)));
elseif index==7
% norm plot
% figure
LW = 'linewidth'; lw = 1.6;
y=chebcoeffs2chebvals(a0);
plot(tspan,real(y(:,N+1)),'r-',LW,lw); hold on
plot(tspan,imag(y(:,N+1)),'b-',LW,lw);
plot(tspan,sum(abs(y),2)-abs(y(:,N+1)),'g-',LW,lw);
% xlabel('$t$','interpreter','latex'), ylabel('$\|u(t)\|_{\infty}$','interpreter', 'latex')
% title('The maximum norm of the solution.')
elseif index==8
% norm plot
% figure
LW = 'linewidth'; lw = 1.6;
y=chebcoeffs2chebvals(a0).'-phi0;
semilogy(tspan,abs(max(y)),LW,lw);
xlabel('$t$','interpreter','latex'), ylabel('$\|u(t)\|_{\infty}$','interpreter', 'latex')
title('The maximum norm of the solution.')
end