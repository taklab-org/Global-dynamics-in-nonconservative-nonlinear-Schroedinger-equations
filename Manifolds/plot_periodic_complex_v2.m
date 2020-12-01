function [real_u,imag_u,x] = plot_periodic_complex_v2(a)

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

real_u = real(u);
imag_u = imag(u);

end

