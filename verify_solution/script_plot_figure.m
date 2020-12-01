LW = 'linewidth';

% fig5a-1
% load data_GE_NLS.mat
figure
plot(mid([y(:,1),y(:,2)])',mid(abs([y(:,6),y(:,6)]))','Color',[0 0.4470 0.7410],LW,1.6)
% hold on
% plot(mid([y(:,1),y(:,2)])',mid([y(:,6)+y(:,8),y(:,6)+y(:,8)])','r',LW,1.6)
% plot(mid([y(:,1),y(:,2)])',mid([y(:,6)-y(:,8),y(:,6)-y(:,8)])','r',LW,1.6)
% stairs(mid([y(:,1)]),mid(y(:,6)),LW,1.6)
xlabel('$t$','interpreter','latex')
ylabel('$\|\bar{a}\|$','interpreter', 'latex')
axis([0,tspan(2),0,1.05*max(mid(y(:,6)))])
title('The values of $\sup_{t}\|\bar{a}(t)\|_{\ell^1}$ in each time step','interpreter', 'latex')
grid on
% % % fig5a-2

figure
px = real(mid(y(:,7))); py = imag(mid(y(:,7)));
plot(px,py,'r',LW,1.6)
hold on
% quiver(px(1:end-1),py(1:end-1),diff(px),diff(py),LW,1.6)
% plot(px(1),py(1),'ro','MarkerSize',10,...
%     'MarkerEdgeColor','red',...
%     'MarkerFaceColor',[1 .6 .6])
xlabel('$\mathrm{Re}\int_0^1 u(t,x)dx$','interpreter','latex')
ylabel('$\mathrm{Im}\int_0^1 u(t,x)dx$','interpreter', 'latex')
title('Solution behavior along with the homogeneous dynamics','interpreter', 'latex')
grid on

z = sup(norm_phi0)*exp(1i*pi*(0:0.05:1));
rz = real(z);
iz = imag(z);
% plot([rz,fliplr(rz)],[iz,zeros(size(rz))],LW,1.6);
fill([rz,fliplr(rz)],[iz,zeros(size(rz))],[1 .6 .6],'EdgeColor','r')

% figure
% semilogy(mid([y(:,1);y(end,2)]),mid([y(:,8);err_at_endpoint]),'.-',LW,1.6)
% xlabel('$t$','interpreter','latex')
% ylabel('$\varepsilon_i$','interpreter', 'latex')
% yticks([1e-14,1e-12,1e-10,1e-8])
% % xticklabels({
% %   '0','10^{-16}','10^{-14}','10^{-12}','10^{-10}',...
% %   '10^{-8}','10^{-6}','10^{-4}','10^{-2}','1'
% % })
% axis([0,0.21,1e-14,2*max(mid([y(:,8);err_at_endpoint]))])
% title('The point-wise error estimate ($\theta=\pi/3$)','interpreter', 'latex')
% % 
% % fig5b-1
% load data_GE_45.mat
% figure
% plot(mid([y(:,1),y(:,2)])',mid([y(:,6),y(:,6)])','Color',[0 0.4470 0.7410],LW,1.6)
% % stairs(mid([y(:,1)]),mid(y(:,6)),LW,1.6)
% xlabel('$t$','interpreter','latex')
% ylabel('$\|\bar{a}\|_X$','interpreter', 'latex')
% axis([0,0.155,0,1.05*max(mid(y(:,6)))])
% title('The values of $\bar{a}$ in $X$-norm ($\theta=\pi/4$)','interpreter', 'latex')
% % fig5b-2
% figure
% semilogy(mid([y(:,1);y(end,2)]),mid([y(:,8);err_at_endpoint]),'.-',LW,1.6)
% xlabel('$t$','interpreter','latex')
% ylabel('$\varepsilon_i$','interpreter', 'latex')
% yticks([1e-14,1e-12,1e-10,1e-8,1e-6])
% axis([0,0.155,1e-14,2*max(mid([y(:,8);err_at_endpoint]))])
% title('The point-wise error estimate ($\theta=\pi/4$)','interpreter', 'latex')
% 
% % fig5c-1
% load data_GE_30.mat
% figure
% plot(mid([y(:,1),y(:,2)])',mid([y(:,6),y(:,6)])','Color',[0 0.4470 0.7410],LW,1.6)
% xlabel('$t$','interpreter','latex')
% ylabel('$\|\bar{a}\|_X$','interpreter', 'latex')
% axis([0,0.132,0,1.05*max(mid(y(:,6)))])
% title('The values of $\bar{a}$ in $X$-norm ($\theta=\pi/6$)','interpreter', 'latex')
% % fig5c-2
% figure
% semilogy(mid([y(:,1);y(end,2)]),mid([y(:,8);err_at_endpoint]),'.-',LW,1.6)
% xlabel('$t$','interpreter','latex')
% ylabel('$\varepsilon_i$','interpreter', 'latex')
% yticks([1e-14,1e-12,1e-10,1e-8,1e-6,1e-4])
% axis([0,0.132,1e-14,2*max(mid([y(:,8);err_at_endpoint]))])
% title('The point-wise error estimate ($\theta=\pi/6$)','interpreter', 'latex')
% 
% % fig5d-1
% load data_GE_15.mat
% figure
% plot(mid([y(:,1),y(:,2)])',mid([y(:,6),y(:,6)])','Color',[0 0.4470 0.7410],LW,1.6)
% xlabel('$t$','interpreter','latex')
% ylabel('$\|\bar{a}\|_X$','interpreter', 'latex')
% axis([0,0.125,0,1.05*max(mid(y(:,6)))])
% title('The values of $\bar{a}$ in $X$-norm ($\theta=\pi/12$)','interpreter', 'latex')
% % fig5d-2
% figure
% semilogy(mid([y(:,1);y(end,2)]),mid([y(:,8);err_at_endpoint]),'.-',LW,1.6)
% xlabel('$t$','interpreter','latex')
% ylabel('$\varepsilon_i$','interpreter', 'latex')
% yticks([1e-12,1e-8,1e-4])
% axis([0,0.125,1e-14,2*max(mid([y(:,8);err_at_endpoint]))])
% title('The point-wise error estimate ($\theta=\pi/12$)','interpreter', 'latex')
