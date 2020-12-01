function norm_of_a = norma(a,nu)

m = length(a);

nu_power = [1 2*nu.^abs(1:m-1)]';
norm_of_a = sum(abs(a).*nu_power);

end