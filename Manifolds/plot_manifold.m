function [] = plot_manifold(p,M,N)

p = reshape(p,N+1,M+1);

figure
hold on

sigma = (-1:.1:1);

n = length(sigma);

for j = 1 : n
     P_at_sigma_j = sum(p.*(repmat(sigma(j).^(0:M),N+1,1)),2);% one-sided Fourier
     plot_periodic_complex(P_at_sigma_j)
end

plot_steady_state(p(:,1))

set(gca,'FontSize',20)

axis tight
xlabel('$$x$$', 'Interpreter', 'latex', 'FontSize', 30)
ylabel('$$Re(u),Im(u)$$', 'Interpreter', 'latex', 'FontSize', 30)

end

