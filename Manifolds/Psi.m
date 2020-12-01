function Psi = Psi(a,inu)

m = length(a);
Psi = zeros(m,1);
a_ext = [a;zeros(m-1,1)]; 

for k=1:m-1
    k_prime = (m:k+m)';
    Psi(k+1)=(1/2)*max(mag(abs( a_ext(abs(k-k_prime)+1) )./inu.^abs(k_prime)));
end

end