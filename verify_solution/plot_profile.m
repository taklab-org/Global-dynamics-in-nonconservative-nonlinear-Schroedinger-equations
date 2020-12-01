function plot_profile(a,color)
if length(color)==1
  color(2) = color(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input is the one-sided Chebyshev  in time %%%
%%%                Fourier cosine  in space                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chebfunpref.setDefaults('factory');
% index = 1;
[~,M] = size(a); %M=2*N+1
N_pad = 150;
M = M + 2*N_pad;

N = (M-1)/2;
dx = 1/(2*N-1);
x = dx*(0:2*N);
a0 = [zeros(1,N_pad),a,zeros(1,N_pad)];
% a0 = [a(1,:);2*a(2:end,:)];
a0_max = sum(abs(a0));

% Plot profile:
% figure
subplot(1,2,1);
plot(x,(2*N+1)*real(ifft(ifftshift(chebcoeffs2chebvals(a0),2).')).',color(1),'linewidth',1.6)
% hold on
xlabel('$x$','interpreter','latex'), ylabel('$\mathrm{Re}(\bar{u})$','interpreter', 'latex')%, zlabel('$\mathrm{Re}(u)$','interpreter', 'latex')
% title('Real part')
% figure
axis([0,1,-a0_max,a0_max])

subplot(1,2,2);
plot(x,(2*N+1)*imag(ifft(ifftshift(chebcoeffs2chebvals(a0),2).')).',color(2),'linewidth',1.6)
% hold on
xlabel('$x$','interpreter','latex'), ylabel('$\mathrm{Im}(\bar{u})$','interpreter', 'latex')%, zlabel('$\mathrm{Im}(u)$','interpreter', 'latex')
% title('Imaginary part')
axis([0,1,-a0_max,a0_max])
