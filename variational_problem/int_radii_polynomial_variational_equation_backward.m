function [r_minus,success] = int_radii_polynomial_variational_equation_backward(C,ta,h,N,n,m,n_pad_cheb,nu)

success = 0 ;

%[Fourier_tail_ta,Cheb_tail_ta] = eval_tails(reshape(ta,n,2*N+1));

n_data = n;
ta0 = ta;
ta = padding_cheb(ta,n,N,n_pad_cheb);
n = n_data + n_pad_cheb;

inu = intval(nu); 

% disp(['nu1 = ',num2str(nu(1)),', n_data = ',num2str(n_data),', N = ',num2str(m),', n_proof = ',num2str(n)])

% disp(['total size = ',num2str(n*(2*m+1))])

% disp('Computing derivative...')
ita = intval(ta); %iangle = intval(angle); 
ih = intval(h);
iDF = int_DF_Fourier_Chebyshev_backward(ita,ih,N,n,m);
DF = mid(iDF);

% disp('Computing inverse...')
A = DF\eye(n*(2*m+1));
iA = intval(A);

% disp('Computing Y0...')
Y0 = zeros(2*m+1,1); % there is a Y0 for each 2N+1 IVP

I = eye(2*m+1);

for j = 1:2*m+1
    
    b = I(:,j);
    c = C(:,j);
    c = padding_cheb(c,n_data,m,n_pad_cheb);
    c = reshape(c,n,2*m+1);
    
    F_ext = int_F_Fourier_Chebyshev_backward_ext(c,ita,b,h,N,n,m);
    F_ext = reshape(F_ext,2*n-1,2*m+1);
    F_finite = F_ext(1:n,:); F_finite = reshape(F_finite,n*(2*m+1),1);
    AF_finite = iA*F_finite; AF_finite = reshape(AF_finite,n,2*m+1);
    AF_tail = diag(1./intval(n:2*n-2))*F_ext(n+1:end,:);
    AF = [AF_finite;AF_tail];
    Y0(j) = sup(compute_nu_norm(AF,inu));
end

Y0 = max(Y0);

% disp(['Y0 = ',num2str(Y0)])

% disp('Computing Z0...')

B = intval(eye(n*(2*m+1))) - iA*iDF;
Z0 = sup(compute_Z0(B,m,n,inu));
% disp(['Z0 = ',num2str(Z0)])

% disp('Computing Z1...')

ih = intval(h); ipi = intval('pi'); in = intval(n);

Psi = zeros(n,2*m+1);
ta_ext = [zeros(n,N) reshape(ita,n,2*N+1) zeros(n,N)];

for k=-m:m
    for ell = 1:n-1        
        k2 = (-m:m); ell_2 = (n:ell+n-1)';
        weights = (inu(1).^abs(ell_2))*(inu(2).^abs(k2));
        Psi(ell+1,k+m+1) = max(max(sup(abs(ta_ext(ell_2-ell,k-k2+2*N+1))./weights)));  
    end
end

T = diag(ones(1,n-1),-1) + diag(ones(1,n-1),1) ; T(1,2)=0;

z_F = h*abs(T)*intval(Psi); z_F(1,:) = 2/nu(1)^n; z_F = reshape(z_F,n*(2*m+1),1);
v_F = abs(iA)*z_F; 
v_F = reshape(v_F,n,2*m+1);
lambda_m = abs(-(ih/2)*(2*ipi)^2*m^2);
z1_finite = compute_nu_norm(v_F,inu);
ta0 = reshape(ta0,n_data,2*N+1);
z1_tail = (1/(2*in))*(inu(1)+1/inu(1))*(abs(lambda_m)+4*ih*compute_nu_norm(ta0,inu));
%Z1 = max(z1_finite,z1_tail);
Z1 = sup(z1_finite+z1_tail);
% disp(['z1_finite = ',num2str(sup(z1_finite)),', z1_tail = ',num2str(sup(z1_tail))])

if sup(intval(Z0)+intval(Z1))<1
    success=1;
end

r_minus = sup(intval(Y0)/(1-intval(Z0)-intval(Z1)));

end

