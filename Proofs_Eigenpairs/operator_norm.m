function op_norm = operator_norm(C,inu,tail)

m = (length(C)-1)/2;
i0 = 1; i1 = (2:m+1); i2 = (m+2:2*m+1);

C00 = C(i0,i0); C01 = C(i0,i1); C02 = C(i0,i2);
C10 = C(i1,i0); C11 = C(i1,i1); C12 = C(i1,i2);
C20 = C(i2,i0); C21 = C(i2,i1); C22 = C(i2,i2);

norm_C00 = abs(C00);
norm_C01 = dual_norm(C01,inu);
norm_C02 = dual_norm(C02,inu); 

norm_C10 = norma(C10,inu);
norm_C20 = norma(C20,inu);

omega = 2*intval('pi');
mu_m = -m.^2*omega^2;

norm_C11 = max([sup(mnorma(m-1,C11,inu)) sup(isequal(tail,1)/abs(mu_m))]);
norm_C12 = mnorma(m-1,C12,inu);
norm_C21 = mnorma(m-1,C21,inu);
norm_C22 = max([sup(mnorma(m-1,C22,inu)) sup(isequal(tail,1)/abs(mu_m))]);

n00 = norm_C00 + norm_C01 + norm_C02; 
n01 = norm_C10 + norm_C11 + norm_C12; 
n02 = norm_C20 + norm_C21 + norm_C22;

op_norm = intval(max([mag(n00) mag(n01) mag(n02)]));

end

