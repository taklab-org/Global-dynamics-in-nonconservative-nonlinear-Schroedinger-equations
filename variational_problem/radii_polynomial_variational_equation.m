function [r_minus,success] = radii_polynomial_variational_equation(C,ta,h,N,n,n_pad_cheb,nu,core)

success = 0 ;
nu1 = nu(1); nu2 = nu(2);

%[Fourier_tail_ta,Cheb_tail_ta] = eval_tails(reshape(ta,n,2*N+1));

n_data = n;
ta0 = ta;
ta = padding_cheb(ta,n,N,n_pad_cheb);
n = n_data + n_pad_cheb;

disp(['nu1 = ',num2str(nu1),', n_data = ',num2str(n_data),', N = ',num2str(core),', n_proof = ',num2str(n)])

disp(['total size = ',num2str(n*(2*core+1))])

disp('Computing derivative...')
DF = DF_Fourier_Chebyshev(ta,h,N,n,core);% interval

disp('Computing inverse...')
A = DF\eye(n*(2*core+1));

disp('Computing Z0...')

% The computation of Z0 still needs to be fixed by incorprating nu1>1
B = eye(n*(2*core+1)) - A*DF;% interval
Z0 = 2*norm(B,1);

disp('Computing Y0...')
Y0 = zeros(2*core+1,1); % there is a Y0 for each 2N+1 IVP

%nu_var_equation = 2;
% weight = nu_var_equation.^-abs(-N:N)';

B = eye(2*core+1);

for j = 1:2*core+1
    
    b = B(j,:);
%     b = weight(j)*B(j,:);
    c = C(:,j);
    %[Fourier_tail,Cheb_tail] = eval_tails(reshape(c,n_data,2*N+1))
    c = padding_cheb(c,n_data,core,n_pad_cheb);
    c = reshape(c,n,2*core+1);
    F_ext = F_Fourier_Chebyshev_ext(c,ta,b,h,N,n,core);% interval
    F_ext = reshape(F_ext,2*n-1,2*core+1);
    F_finite = F_ext(1:n,:); F_finite = reshape(F_finite,n*(2*core+1),1);
    AF_finite = A*F_finite; AF_finite = reshape(AF_finite,n,2*core+1);% interval
    AF_tail = diag(1./(n:2*n-2))*F_ext(n+1:end,:);% interval
    AF = [AF_finite;AF_tail];
    Y0(j) = compute_nu_norm(AF,nu);% interval
end

Y0 = max(Y0);

disp(['Y0 = ',num2str(Y0)])

disp('Computing Z1...')


Psi = zeros(n,2*core+1);
ta_ext = [zeros(n,N) reshape(ta,n,2*N+1) zeros(n,N)];

% N = core;

for k=-core:core
    for ell = 1:n-1        
        k2 = (-core:core); ell_2 = (n:ell+n-1)';
        weights = (nu1.^abs(ell_2))*(nu2.^abs(k2));% interval
        Psi(ell+1,k+core+1) = max(max(abs(ta_ext(ell_2-ell,k-k2+2*core+1))./weights));  
    end
end

T = diag(ones(1,n-1),-1) + diag(ones(1,n-1),1) ; T(1,2)=0;

z_F = h*T*Psi; z_F(1,:) = 2/nu1^n; z_F = reshape(z_F,n*(2*core+1),1);% interval
v_F = abs(A)*z_F; 
v_F = reshape(v_F,n,2*core+1);
lambda_N = abs(-(h/2)*(2*pi)^2*core^2);% interval
z1_finite = compute_nu_norm(v_F,nu);% interval
ta0 = reshape(ta0,n_data,2*N+1);
z1_tail = (1/n)*(abs(lambda_N)+h*compute_nu_norm(ta0,nu));% interval
%Z1 = max(z1_finite,z1_tail);
Z1 = z1_finite+z1_tail;% interval
disp(['z1_finite = ',num2str(z1_finite),', z1_tail = ',num2str(z1_tail)])

if Z0+Z1<1
    success=1;
end

r_minus = Y0/(1-Z0-Z1);% interval

end

