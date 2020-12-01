function plot_manifold_v2(p,M,N,lambda,sigma,conj)
arguments
  p;M;N;lambda;sigma;conj=false;
end

% sigma is in {-1,1} (determines the side)
p = reshape(p,N+1,M+1);

% number of accuracy digit in parameter space
N_digits = 13; %3:P_at_1, 4:pt2, 
t_min = -N_digits*log(10)/real(lambda);

t_mesh = (t_min:abs(t_min)/500:0); 
x = (0:.001:1);

n1 = length(t_mesh); n2 = length(x);

surf_real_manifold = zeros(n1,n2);
surf_imag_manifold = zeros(n1,n2);

for j = 1 : n1
    
    t_sigma = exp(lambda*t_mesh(j))*sigma;
    P_at_t_sigma = sum(p.*(repmat(t_sigma.^(0:M),N+1,1)),2);
    [real_u,imag_u] = plot_periodic_complex_v2(P_at_t_sigma);
    surf_real_manifold(j,:) = real_u;
    surf_imag_manifold(j,:) = imag_u;

end

if conj
  subplot(2,1,1)
  % f1 = figure;
  surf(x,-t_mesh,surf_real_manifold,'EdgeColor','none'), hold on
  surf(x+1,-t_mesh,surf_real_manifold,'EdgeColor','none');
  surf(x-1,-t_mesh,surf_real_manifold,'EdgeColor','none');
  surf(x-2,-t_mesh,surf_real_manifold,'EdgeColor','none');
  set(gca,'FontSize',20)
  axis tight
  xlabel('$$x$$', 'Interpreter', 'latex', 'FontSize', 30)
  ylabel('$$t$$', 'Interpreter', 'latex', 'FontSize', 30)
  % zlabel('$$Re(u)$$', 'Interpreter', 'latex', 'FontSize', 30)
  
  subplot(2,1,2)
  % f2 = figure;
  surf(x,-t_mesh,-surf_imag_manifold,'EdgeColor','none'), hold on
  surf(x+1,-t_mesh,-surf_imag_manifold,'EdgeColor','none')
  surf(x-1,-t_mesh,-surf_imag_manifold,'EdgeColor','none')
  surf(x-2,-t_mesh,-surf_imag_manifold,'EdgeColor','none')
  set(gca,'FontSize',20)
  axis tight
  xlabel('$$x$$', 'Interpreter', 'latex', 'FontSize', 30)
  ylabel('$$t$$', 'Interpreter', 'latex', 'FontSize', 30)
  % zlabel('$$Im(u)$$', 'Interpreter', 'latex', 'FontSize', 30)
  return
end

subplot(2,1,1)
% f1 = figure;
surf(x,t_mesh,surf_real_manifold,'EdgeColor','none'), hold on
surf(x+1,t_mesh,surf_real_manifold,'EdgeColor','none');
surf(x-1,t_mesh,surf_real_manifold,'EdgeColor','none');
surf(x-2,t_mesh,surf_real_manifold,'EdgeColor','none');
set(gca,'FontSize',20)
axis tight
xlabel('$$x$$', 'Interpreter', 'latex', 'FontSize', 30)
ylabel('$$t$$', 'Interpreter', 'latex', 'FontSize', 30)
% zlabel('$$Re(u)$$', 'Interpreter', 'latex', 'FontSize', 30)

subplot(2,1,2)
% f2 = figure;
surf(x,t_mesh,surf_imag_manifold,'EdgeColor','none'), hold on
surf(x+1,t_mesh,surf_imag_manifold,'EdgeColor','none')
surf(x-1,t_mesh,surf_imag_manifold,'EdgeColor','none')
surf(x-2,t_mesh,surf_imag_manifold,'EdgeColor','none')
set(gca,'FontSize',20)
axis tight
xlabel('$$x$$', 'Interpreter', 'latex', 'FontSize', 30)
ylabel('$$t$$', 'Interpreter', 'latex', 'FontSize', 30)
% zlabel('$$Im(u)$$', 'Interpreter', 'latex', 'FontSize', 30)

end