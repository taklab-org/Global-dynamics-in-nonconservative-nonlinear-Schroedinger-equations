function x = newton(x,v,theta)

tol = 5e-15; %% tolerance for Newton's method

F = F_eigenpairs(x,v,theta);

nF = norm(F,1);
display(['At the beginning ||F|| = ',num2str(nF)])

k=0;

while (k<=60) && (nF > tol)
    DF = DF_eigenpairs(x,v,theta);
    x = x - DF\F;
    F = F_eigenpairs(x,v,theta);
    nF = norm(F);
    display(['||F|| = ',num2str(nF),', ||DF^(-1)|| = ',num2str(norm(inv(DF),1))])
    k = k+1;                            
end

end