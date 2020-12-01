function plot_glue_manifold(a,tspan,f,conj,N_pad,n_pad)
arguments
  a = []; tspan=[];
  f = figure; conj = false;
  N_pad = 100; n_pad = 100;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input is the one-sided in Chebyshev %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chebfunpref.setDefaults('factory');
% index = 1;
[n,M] = size(a); %M=2*N+1
N = (M-1)/2;
% 
a0 = [a(1,:);2*a(2:end,:);zeros(n_pad,M)];
a0(:,1:N) = fliplr(a0(:,N+2:end));
% a0(:,N+2:end) = fliplr(a0(:,1:N)); % symmetry-condition
% 
N = N + N_pad; n = n + n_pad;
dx = 1/(2*N-1);
x = dx*(0:2*N);
% 
a0 = [zeros(n,N_pad),a0,zeros(n,N_pad)];
% 
if length(tspan)<2
  tspan=chebpts(n,[0,tspan]);
else
  tspan=chebpts(n,tspan);
end

% Plot profile:
if conj
  figure(f);
  subplot(2,1,1)
  hold on
  zr = (2*N+1)*real(ifft(ifftshift(chebcoeffs2chebvals(a0),2).')).';
  surf(x,-tspan,zr,'EdgeColor','none')
  surf(x+1,-tspan,zr,'EdgeColor','none')
  surf(x-2,-tspan,zr,'EdgeColor','none')
  surf(x-1,-tspan,zr,'EdgeColor','none')
  % 
  subplot(2,1,2);
  hold on
  zi = (2*N+1)*imag(ifft(ifftshift(chebcoeffs2chebvals(a0),2).')).';
  surf(x,-tspan,-zi,'EdgeColor','none')
  surf(x+1,-tspan,-zi,'EdgeColor','none')
  surf(x-1,-tspan,-zi,'EdgeColor','none')
  surf(x-2,-tspan,-zi,'EdgeColor','none')
  return
end
figure(f);
subplot(2,1,1)
hold on
% subplot(1,2,1);
zr = (2*N+1)*real(ifft(ifftshift(chebcoeffs2chebvals(a0),2).')).';
surf(x,tspan,zr,'EdgeColor','none')
surf(x+1,tspan,zr,'EdgeColor','none')
surf(x-2,tspan,zr,'EdgeColor','none')
surf(x-1,tspan,zr,'EdgeColor','none')
% xlabel('$x$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Re}(u)$','interpreter', 'latex')
% title('Real part')
% figure(f2);
subplot(2,1,2);
hold on
zi = (2*N+1)*imag(ifft(ifftshift(chebcoeffs2chebvals(a0),2).')).';
surf(x,tspan,zi,'EdgeColor','none')
surf(x+1,tspan,zi,'EdgeColor','none')
surf(x-1,tspan,zi,'EdgeColor','none')
surf(x-2,tspan,zi,'EdgeColor','none')
% xlabel('$x$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Im}(u)$','interpreter', 'latex')
% title('Imaginary part')