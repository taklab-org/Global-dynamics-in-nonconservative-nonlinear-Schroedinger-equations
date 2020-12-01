function [I,success] = rad_poly_eigenpairs(x,v,theta,nu)

success = 0; I = [];

DF = DF_eigenpairs(x,v,theta);
A = inv(DF);
iA = intval(A);

m = (length(x)-1)/2;

i0 = 1; i1 = (2:m+1); i2 = (m+2:2*m+1);

lambda = x(i0);
a = x(i1); %%% Fourier coefficients of the steady state
b = x(i2); %%% Fourier coefficients of the eigenfunction 

%%%% Y0 %%%%

a_ext = [a;zeros(m-1,1)]; 
b_ext = [b;zeros(m-1,1)];
v_ext = [v;zeros(m-1,1)];

ix_ext = intval([lambda;a_ext;b_ext]);
inu = intval(nu);

F_ext = iF_eigenpairs(ix_ext,v_ext,theta);
F0_ext = F_ext(1);
F1_ext = F_ext(2:2*m);
F2_ext = F_ext(2*m+1:4*m-1);

omega = 2*intval('pi');
k_tail = (m:2*m-2)';
mu_k_tail = -k_tail.^2*omega^2;

y_F = abs(iA*[F0_ext;F1_ext(1:m);F2_ext(1:m)]);
y0_F = y_F(1);
y1_F = y_F(2:m+1); 
y2_F = y_F(m+2:2*m+1);
y1_tail = abs(F1_ext(m+1:2*m-1)./mu_k_tail); 
y2_tail = abs(F2_ext(m+1:2*m-1)./mu_k_tail);
y0 = y0_F;
y1 = [y1_F;y1_tail]; 
y2 = [y2_F;y2_tail];

nu_power = [1 2*inu.^(1:2*m-2)]';

Y0 = intval(max([mag(y0) mag(sum(y1.*nu_power)) mag(sum(y2.*nu_power))] ));

%%%% Z0 %%%%

iDF = iDF_eigenpairs(x,v,theta);

B = intval(eye(2*m+1)-iA*iDF);

Z0 = operator_norm(B,inu,0);

%%%% Z1 %%%%

hQa = zeros(m,1); hQb = zeros(m,1);

for k=1:m-1
    k_prime = (m:k+m)';
    hQa(k+1)=(1/2)*max(mag(abs( a_ext(abs(k-k_prime)+1) )./inu.^abs(k_prime)));
    hQb(k+1)=(1/2)*max(mag(abs( b_ext(abs(k-k_prime)+1) )./inu.^abs(k_prime)));
end

hQa = intval(hQa);
hQb = intval(hQb);

iA01 = abs(iA(i0,i1)); iA02 = abs(iA(i0,i2)); 
iA11 = abs(iA(i1,i1)); iA12 = abs(iA(i1,i2)); 
iA21 = abs(iA(i2,i1)); iA22 = abs(iA(i2,i2));

imu_m = 4*m^2*intval('pi')^2;

Z1_0 = 2*iA01*hQa + 2*iA02*(hQa+hQb);
Z1_1 = 2*norma(iA11*hQa,inu) + 2*norma(iA12*(hQa+hQb),inu) + norma(a,inu)/imu_m;
Z1_2 = 2*norma(iA21*hQa,inu) + 2*norma(iA22*(hQa+hQb),inu) + (2*norma(a,inu)+2*norma(b,inu)+abs(lambda))/imu_m;

Z1 = intval(max([mag(Z1_0) mag(Z1_1) mag(Z1_2)]));

%%%% Z2 %%%%

normA = operator_norm(iA,inu,1);

Z2 = 6*normA;

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