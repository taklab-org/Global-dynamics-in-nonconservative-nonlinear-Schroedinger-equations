function [b_end,tot_I] = integrate_Fisher(b,lambda,h,bar_j,steps)

bar_k = length(b);

plot_b(b)
hold on

A=zeros(bar_k,bar_j);
a=reshape(A,bar_k*bar_j,1);

tol=1e-14;

F = F_Fourier_Chebyshev_Fisher(a,b,lambda,h,bar_k,bar_j);
display(['Before step #1, ||F||_1 = ',num2str(norm(F,1))])
itercount=0;

L = linear_operator_L_fisher(h,bar_k,bar_j);
L_inv = inv(L);
norm_inv_L=norm(L_inv,1);
display(['||inv(DF)||_1 = ',num2str(norm_inv_L)])

tot_I = zeros(steps,2);

for k=1:steps
    
    disp(' ------------------------' )
    disp(['  Integration step = ',num2str(k)])
    disp(' ------------------------' )
    
    F = F_Fourier_Chebyshev_Fisher(a,b,lambda,h,bar_k,bar_j);
    norm_F = norm(F,1);
    
    while (itercount<=100) && (norm_F > tol),
        
        a=a-L_inv*F;
        F = F_Fourier_Chebyshev_Fisher(a,b,lambda,h,bar_k,bar_j);
        norm_F = norm(F,1);
        display(['At Newton iteration # ',num2str(itercount+1),', ||F||_1 = ',num2str(norm(F,1))])
        itercount=itercount+1;
        
    end
    
    tail_norms(a,bar_k,bar_j);

    nu=1.01;
    
    [I] = radii_polynomial_integrate_Fisher(a,b,lambda,bar_k,bar_j,h,nu);
    
    tot_I(k,:) = I;
    
    itercount=0;
    
    A = reshape(a,bar_j,bar_k)';
    a_end = zeros(bar_k,1);
    
    for index=1:bar_k
        a_end(index,1)= A(index,1)+2*sum(A(index,2:end));
    end
    
    b=b+a_end;
    
    A=zeros(bar_k,bar_j);
    a=reshape(A,bar_k*bar_j,1);
    
end

b_end = b;

plot_b(b_end)

end

