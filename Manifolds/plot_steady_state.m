function [] = plot_steady_state(a)

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

plot(x,real(u),'black','linewidth',3)
plot(x,imag(u),'black','linewidth',3)

end

