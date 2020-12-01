function [Df] = finite_diff_DF(p,par)

h = 1e-6;
m = length(p);
E = eye(m);
Df = zeros(m);

for j = 1:m
    ph = p + h*E(:,j);
    Df(:,j) = (f_1d_unstable_manifold(ph,par)-f_1d_unstable_manifold(p,par))/h;
end

end