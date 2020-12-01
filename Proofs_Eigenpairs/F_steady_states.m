function F = F_steady_states(a,theta)

N = length(a)-1;
omega = 2*pi;
k = (0:N)';

a_ext = [flip(a(2:end),1);a];

a2 = quadratic(a_ext,a_ext); 
F = exp(1i*theta)*(-k.^2*omega^2.*a + a2(N+1:end)) ; 

end

