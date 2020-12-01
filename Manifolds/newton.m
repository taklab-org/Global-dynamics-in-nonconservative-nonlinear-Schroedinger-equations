function p = newton(p,par)

tol = 5e-11; %% tolerance for Newton's method

f = f_1d_unstable_manifold(p,par);

nf = norm(f,1);
display(['At the beginning ||f|| = ',num2str(nf)])

k=0;

while (k<=10) && (nf > tol)
    Df = Df_1d_unstable_manifold(p,par);
    p = p - Df\f;
    f = f_1d_unstable_manifold(p,par);
    nf = norm(f);
    display(['||f|| = ',num2str(nf),', ||Df^(-1)|| = ',num2str(norm(inv(Df),1))])
    k = k+1;                            
end

end