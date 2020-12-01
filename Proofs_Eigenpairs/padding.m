function [x,v] = padding(x,v,n)

m = (length(x)-1)/2;

a = [x(2:m+1);zeros(n,1)];
b = [x(m+2:2*m+1);zeros(n,1)];

x = [x(1);a;b];
v = [v;zeros(n,1)];

end

