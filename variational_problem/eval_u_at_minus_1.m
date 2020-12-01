function u_at_minus_1 = eval_u_at_minus_1(c)

n = length(c);

u_at_minus_1 = c(1) + 2*sum( ((-1).^(1:n-1)').*c(2:n)) ; 

end

