function [I,success] = rad_poly_1d_unstable_CGL(p,par,r0,nu)

success = 0; I = [];

Df = Df_1d_unstable_manifold(p,par);
A = inv(Df);
iA = intval(A);

nu = intval(nu);

%%%%%%%%%%%%
%%%% Y0 %%%%
%%%%%%%%%%%%

N = par(1); M = par(2);
f_ext = int_f_1d_unstable_manifold_ext(p,par,r0);
f_ext = reshape(f_ext,2*N+1,2*M+1);
f_F = f_ext(1:N+1,1:M+1); f_F = reshape(f_F,(N+1)*(M+1),1);

omega = 2*intval('pi');

%Y_0^(1)%

matrix_m = repmat((0:2*M),2*N+1,1);
matrix_k = repmat((0:2*N)',1,2*M+1);

theta = intval(par(2*N+6));
lambda = midrad(par(N+4),r0);
mu_km = lambda*matrix_m + exp(1i*theta)*matrix_k.^2*omega^2;

k = abs(0:2*N)'; weights = (nu.^k)*ones(1,2*M+1);

y01_1 = sum(sum(abs(f_ext(N+2:2*N+1,3:M+1)./mu_km(N+2:2*N+1,3:M+1)).*weights(N+2:2*N+1,3:M+1)));
y01_2 = sum(sum(abs(f_ext(1:2*N+1,M+2:2*M+1)./mu_km(1:2*N+1,M+2:2*M+1)).*weights(1:2*N+1,M+2:2*M+1)));

Y0_1 = compute_nu_norm(iA*f_F,par,nu) ...
    + norma([zeros(N+1,1);f_ext(N+2:2*N+1,1)],nu) ...
    + norma([zeros(N+1,1);f_ext(N+2:2*N+1,2)],nu) ...
    + y01_1 + y01_2;

sup(Y0_1)

%Y_0^(2)%

k = (0:N);

alpha_0 = 0; alpha_1 = 0;
for j = 0:M
    A_j0 = A(k+1+j*(N+1),k+1);
    A_j1 = A(k+1+j*(N+1),k+1+(N+1));
    alpha_0 = alpha_0 + mnorma(A_j0,N,isequal(j,0),nu);
    alpha_1 = alpha_1 + mnorma(A_j1,N,isequal(j,0),nu);
end

Y0_2 = (alpha_0 + alpha_1)*r0;

Y0 = sup(Y0_1 + Y0_2);
disp(['Y0 = ',num2str(Y0)])

%%%%%%%%%%%%
%%%% Z0 %%%%
%%%%%%%%%%%%

disp('Computing Z0...')

iDf = iDf_1d_unstable_manifold(p,par,r0);
B = intval(eye((N+1)*(M+1))) - iA*iDf;
Z0 = sup(operator_norm(B,par,0,nu));

disp(['Z0 = ',num2str(Z0)])

%%%%%%%%%%%%
%%%% Z1 %%%%
%%%%%%%%%%%%

p = reshape(p,N+1,M+1);

PSI = zeros(N+1,M+1);

for ell = 1:M+1
    PSI(:,ell) = Psi(p(:,ell),nu);
end

PSI = intval(PSI);

z1 = intval(zeros(N+1,M+1));

for m = 0:M
    z1(:,m+1) = 2*sum(PSI(:,1:m+1),2);
end

z1(:,1:2)=0;

z1 = reshape(z1,(N+1)*(M+1),1);

Z1_1 = compute_nu_norm(abs(iA)*z1,par,nu);

%%%% Bound for CGL
%Z1_2 = 2*(1/(real(lambda)*(M+1))+1/(cos(theta)*(N+1)^2*omega^2))*compute_nu_norm(p,par,nu);

%%%% Bound for NLS
lbda_r = real(lambda); lbda_i = imag(lambda);
m = (0:M);
if min(inf(lbda_i*m + (N+1)^2*omega^2))<0
    disp('error with the tail bound in mu')
    return
else    
    mu_star = intval(min(inf(sqrt((lbda_r*m).^2 + (lbda_i*m + (N+1)^2*omega^2 ).^2))));
end
    
Z1_2 = 2*( 1/mu_star + 1/(lbda_r*(M+1)) )*compute_nu_norm(p,par,nu);

Z1 = Z1_1 + Z1_2;

disp(['Z1 = ',num2str(sup(Z1))])

%%%%%%%%%%%%
%%%% Z2 %%%%
%%%%%%%%%%%%

normA = operator_norm(iA,par,1,nu);

Z2 = 2*normA;

disp(['Z2 = ',num2str(sup(Z2))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Radii Polynomial Verification %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display(['Y0 = ',num2str(mid(Y0)),', Z0 = ',num2str(mid(Z0)),', Z1 = ',num2str(mid(Z1)),', Z2 = ',num2str(mid(Z2))])

if inf(1-Z0-Z1)>0
    if inf((1-Z0-Z1)^2-4*Y0*Z2) > 0
        rmin=sup(((1-Z0-Z1) - sqrt((1-Z0-Z1)^2-4*Y0*Z2))/(2*Z2));
        rmax=inf(((1-Z0-Z1) + sqrt((1-Z0-Z1)^2-4*Y0*Z2))/(2*Z2));
        if rmin<rmax
            success=1;
            I=[rmin rmax];
        else
            disp('failure: rmin > rmax')
        end
    else
        disp('failure: discriminant is negative')
    end
else
    disp('failure: linear term is positive')
end

end