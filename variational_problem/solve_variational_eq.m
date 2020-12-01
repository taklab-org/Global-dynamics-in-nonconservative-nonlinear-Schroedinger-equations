function C = solve_variational_eq(ta,angle,h,N,n)

tol = 5e-14;

DF = DF_Fourier_Chebyshev(ta,angle,h,N,n);

DF_inv = DF\eye(n*(2*N+1));

delta = 1e-3;

I = delta*eye(2*N+1);
C = zeros(n*(2*N+1));

for j=1:2*N+1
    
    b = I(:,j);
    c = zeros(n,2*N+1); c(1,:) = b;
    c = reshape(c,(2*N+1)*n,1); 
        
    F = F_Fourier_Chebyshev(c,ta,b,angle,h,N,n);
    disp(j)
    display(['Before step #1, ||F||_1 = ',num2str(norm(F,1))])
    itercount=0;
    
    while (itercount<=10) && (norm(F,1) > tol)
        
        c = c - DF_inv*F;
        F = F_Fourier_Chebyshev(c,ta,b,angle,h,N,n);
        display(['After step # ',num2str(itercount+1),', ||F||_1 = ',num2str(norm(F,1))])
        itercount=itercount+1;
        
    end
    
    C(:,j)=c;
    
end

end
