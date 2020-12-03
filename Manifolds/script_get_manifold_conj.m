clear
clc
close all

load eigenpair_pt1_conj_NLS_proof % conjugate of pt1
% load eigenpair_pt2_conj_NLS_proof % conjugate of pt2

% M = Taylor projection;  n_ext = Fourier padding
% conjugate of pt1
M = 150; scaling = 20; nu = 1; n_ext = 10;  % gives end_points_proof2_NLS
% conjugate of pt2
% angle = -60/180*pi;
% M = 60; scaling = 75*exp(1i*angle); nu = 1; n_ext = 10;


N = length(a) - 1 ; % Fourier projection

N = N + n_ext; a = [a;zeros(n_ext,1)]; b = [b;zeros(n_ext,1)];

p = zeros(N+1,M+1);

p(:,1) = a; % The fixed point
p(:,2) = scaling*b; % The eigenvector

par = zeros(2*N+6,1);

par(1) = N; par(2) = M;
par(3:N+3) = a;
par(N+4) = lambda;
par(N+5:2*N+5) = p(:,2);
par(2*N+6) = theta;

p = reshape(p,(N+1)*(M+1),1); 

p = newton(p,par); %p = reshape(p,N+1,M+1);

plot_manifold(p,M,N)

%%% We impose that the first oder data are exact
p = reshape(p,N+1,M+1);
p(:,1) = a;
p(:,2) = scaling*b;

p = reshape(p,(N+1)*(M+1),1);
p = clean_p(p);
[I,success] = rad_poly_1d_unstable_CGL(p,par,r0,nu);
disp(['Success = ',num2str(success)])

rmin = I(1);
p = reshape(p,N+1,M+1);
P_at_1 = sum(intval(p),2); error_at_1 = rmin + norma(rad(P_at_1),nu);
%plot_periodic_complex(P_at_1)
%hold on
P_at_minus_1 = sum(intval(p).*(repmat((-1).^(0:M),N+1,1)),2); 
error_at_minus_1 = rmin + norma(rad(P_at_minus_1),nu);
%plot_periodic_complex(P_at_minus_1)

disp(['error_at_1 = ',num2str(error_at_1),', error_at_minus_1 = ',num2str(error_at_minus_1)])

% save end_points_proof_NLS_pt2 P_at_1 P_at_minus_1 error_at_1 error_at_minus_1
% end