close all
clear

% load eigenpair1
% load eigenpair2
% load eigenpair_NLS
% load eigenpair_NLS_v2 % pt1
% load eigenpair_NLS_pt2 % pt2
% load eigenpair_NLS_pt1_conj % conjugate of pt1
load eigenpair_NLS_pt2_conj % conjugate of pt2

v = clean_p(v); x = clean_p(x); 

nu = 1;

[I,success] = rad_poly_eigenpairs(x,v,theta,nu);

disp('   ')
disp(['I = [',num2str(I(1)),',',num2str(I(2)),']'])

m = (length(x)-1)/2; a = x(2:m+1); plot_periodic_complex(a)

lambda = x(1);
b = x(m+2:2*m+1);
r0 = I(1);

save ../Manifolds/eigenpair_pt2_conj_NLS_proof a b lambda theta r0 nu