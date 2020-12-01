function [] = plot_periodic_complex(a)

N = (length(a)-1);
omega = 2*pi;

x = (0:.001:1);
n = length(x);
k = (-N:N)';

a = [flipud(a(2:end)); a];

u = zeros(1,n);

for j = 1:n
    u(j) = sum(a.*exp(1i*k*omega*x(j)));
end

plot(x,real(u),'blue','linewidth',3)
hold on
plot(x,imag(u),'red','linewidth',3)
hold off

set(gca,'FontSize',20)

axis tight
xlabel('$$x$$', 'Interpreter', 'latex', 'FontSize', 30)
ylabel('$$Re(u),Im(u)$$', 'Interpreter', 'latex', 'FontSize', 30)


end

