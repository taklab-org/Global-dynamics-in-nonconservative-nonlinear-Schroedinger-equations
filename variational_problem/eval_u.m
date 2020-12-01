function [u]=eval_u(x,m1,m2,t,y)

%%% m1 : size of the time Fourier expansion
%%% m2 : size of the spatial Fourier expansion

n=2*m1*m2-2*m1-m2+2;
L=x(1);
A=x(2:m1*m2-m1-m2+2);
B=x(m1*m2-m1-m2+3:n);

a=zeros(m1-1,m2-1);
b=zeros(m1,m2-1);

for k2=1:m2-1
    a(:,k2)=A((k2-1)*(m1-1)+1:k2*(m1-1));
    b(:,k2)=B((k2-1)*m1+1:k2*m1);
end

c=a+1i*b(2:m1,1:m2-1);
    
PSI_1=exp(1i*L*(1:m1-1)'*t)*exp(1i*(1:m2-1)*y);
PSI_2=exp(-1i*L*(1:m1-1)'*t)*exp(1i*(1:m2-1)*y);
PSI_3=exp(-1i*L*(1:m1-1)'*t)*exp(-1i*(1:m2-1)*y);
PSI_4=exp(1i*L*(1:m1-1)'*t)*exp(-1i*(1:m2-1)*y);

b_term=-2*sum(b(1,1:m2-1).*sin((1:m2-1)*y));

u=b_term+real(sum(sum(c.*PSI_1-conj(c).*PSI_2+conj(c).*PSI_3-c.*PSI_4)));

end