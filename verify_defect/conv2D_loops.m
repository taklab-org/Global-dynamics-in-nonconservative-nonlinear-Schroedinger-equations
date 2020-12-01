function [ab] = conv2D_loops(a,b)

% Inputs : a = (a_{ell,k}) and b = (b_{ell,k})
%           ell = 0,...,n (Chebyshev indices)
%           k = -N,...,N (Fourier indices)

[n1,n2] = size(a);

n = n1-1;
N = (n2-1)/2;

ab = zeros(n1,n2);

for ell = 0:n
    for k = -N:N        
        sum = 0 ;       
        for k1 = -N:N
            for ell_1 = -n:n
                if abs(k-k1)<=N && abs(ell-ell_1)<=n
                    sum = sum + a(abs(ell_1)+1,k1+N+1)*b(abs(ell-ell_1)+1,k-k1+N+1);
                end
            end
        end       
        ab(ell+1,k+N+1) = sum;       
    end
end


end

