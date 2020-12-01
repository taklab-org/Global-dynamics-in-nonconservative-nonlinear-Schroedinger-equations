function [P_at_minus_1,P_at_1] = produce_unstable_data_from_angle(angle,extra_dim)

% Input: angle for the eigenvector

% load eigenpair_NLS_proof
load eigenpair_pt2_NLS_proof
a = [a;zeros(extra_dim,1)]; b = [b;zeros(extra_dim,1)];
% M = 50; scaling = 20*exp(1i*angle);
M = 40; scaling = 55*exp(1i*angle);
 M = 150; scaling = 20; nu = 1; n_ext = 10;  % gives end_points_proof2_NLS


N = length(a) - 1 ; % Fourier projection

p = zeros(N+1,M+1);

p(:,1) = a; % The fixed point
p(:,2) = scaling*b; % The eigenvector

% Verify that we have a correct eigenfunction
b_test = p(:,2);
ab = quadraticFFT(a,b_test);
m = length(b_test);
k = (0:m-1)';
omega = 2*pi;
mu = -k.^2*omega^2;
test_eigs = norm(exp(1i*theta)*(mu.*b_test + 2*ab) - lambda*b_test);

if test_eigs>1e-8
    disp('problem with first order data')
    return
end

par = zeros(2*N+6,1);

par(1) = N; par(2) = M;
par(3:N+3) = a;
par(N+4) = lambda;
par(N+5:2*N+5) = p(:,2);
par(2*N+6) = theta;

p = reshape(p,(N+1)*(M+1),1);

p = newton(p,par); %p = reshape(p,N+1,M+1);

% plot_manifold(p,M,N)

%%% We impose that the first order data are exact
p = reshape(p,N+1,M+1);
p(:,1) = a;
p(:,2) = scaling*b;

p = reshape(p,(N+1)*(M+1),1);
p = clean_p(p);

p = reshape(p,N+1,M+1);

P_at_1 = sum(p,2);
P_at_minus_1 = sum(p.*(repmat((-1).^(0:M),N+1,1)),2);

end
