function [ab,abfull] = convtensor(a,b)
% computes the convolution of general Fourier series A and B
% A and B should be tensors of equal size
% works in all dimensions
% does not assume any symmetry
% the second, optional output is the "extended product"

sz = size(a);
dim=length(sz);

N = (sz-1)/2;
M = 2.^nextpow2(4*N+1);

a1 = altzeros(M(1:dim),a(1));
b1 = altzeros(M(1:dim),b(1));
M2=M/2;
M2(M==1)=0; % take care of dimensions with just one element (N=0)
s1 = M2 + 1 - N;
s2 = M2 + 1 + N;
s1f = M2 + 1 - 2*N;
s2f = M2 + 1 + 2*N;
for j=1:dim
  s{j}=s1(j):s2(j);
  sf{j}=s1f(j):s2f(j);
  fl{j}=[M2(j)+1:M(j),1:M2(j)];
end

a1(s{1:dim})=a;
b1(s{1:dim})=b;

a2 = a1(fl{1:dim}); %fftshift
b2 = b1(fl{1:dim}); %fftshift

u = altfftn(a2);
v = altfftn(b2);
uv = reshape(u(:).*v(:),size(u));   
ab2 = altifftn(uv);

ab1 = ab2(fl{1:dim}); %fftshift

ab = ab1(s{1:dim});
abfull = ab1(sf{1:dim});

return



