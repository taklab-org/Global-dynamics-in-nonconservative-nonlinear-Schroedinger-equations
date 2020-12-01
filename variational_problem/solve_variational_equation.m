function [C,r_minus] = solve_variational_equation(a,h,N,n,m,rigorous)

% clear
% load approximate_solution_v2% nu = [1.02;1]; n_pad = 300; N=20; n=150; %(These choices work!)
% load approximate_solution_v3 %nu = [1.02;1]; n_pad = 300; N=15; n=110;(These choices work?)

addpath('../verify_defect/')

ta = reshape(a,(2*N+1)*n,1);

tol = 5e-14;

DF = DF_Fourier_Chebyshev(ta,h,N,n,m);

DF_inv = DF\eye(n*(2*m+1));

B = eye(2*m+1);
C = zeros(n*(2*m+1),2*m+1);

for j=1:2*m+1
    
    b = B(j,:); 
    c = zeros(n,2*m+1); c(1,:) = b;
    c = reshape(c,(2*m+1)*n,1); 
        
    F = F_Fourier_Chebyshev(c,ta,b,h,N,n,m);
    itercount=0;
    
    while (itercount<=10) && (norm(F,1) > tol)
        
        c = c - DF_inv*F;
        F = F_Fourier_Chebyshev(c,ta,b,h,N,n,m);
        itercount=itercount+1;
        
    end
    
    C(:,j)=c;

%     c = reshape(c,n,2*N+1);
%     cc = [c(1,:);2*c(2:end,:)]; % Chebyshev coefficients
%     C2(:,j)=reshape(cc,(2*N+1)*n,1);% contain two-sided Chebyshev coeffs.
% 
%     for k=1:2*N+1
%         [decay_rate,constant] = exp_decay_a_least_square(c(:,k)); rho(k,j)=decay_rate; K(k,j)=constant; 
%     end
    
end

clear B F b c itercount j delta

lambda_m = abs(-(h/2)*(2*pi)^2*m^2);
n_cheb_proof = ceil(lambda_m); n_pad_Z1_tail = n_cheb_proof - n;
% disp(['min total size = ',num2str(n_cheb_proof*(2*m+1))])
% disp(['n_pad_Z1_tail = ',num2str(n_pad_Z1_tail)])

%nu1_least_square = min(min(rho));
%disp(['nu1 least square = ',num2str(nu1_least_square)])
% nu = [1.03;1]; n_pad = 500;  % This one works!
nu = [1.1;1]; n_pad = 50;
if n_pad<n_pad_Z1_tail
    error('the tail estimate will fail!')
end

% z1_tail = (1/(n+n_pad))*(abs(lambda_m)+h*compute_nu_norm(a,nu));
% disp(['z1_tail = ',num2str(z1_tail)])

if rigorous>0
  [r_minus,success] = int_radii_polynomial_variational_equation(C,ta,h,N,n,m,n_pad,nu);
else
  [r_minus,success] = radii_polynomial_variational_equation(C,ta,h,N,n,n_pad,nu,m);
end
% success
if success~=1
  error('variational problem is not verified!')
end