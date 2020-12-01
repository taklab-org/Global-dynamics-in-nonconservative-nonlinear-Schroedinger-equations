function [s]=quadratic(a1,a2)

m=(length(a1)+1)/2;

ta1=[zeros(m,1);a1;zeros(m,1)]; tu1=ifft(ifftshift(ta1));
ta2=[zeros(m,1);a2;zeros(m,1)]; tu2=ifft(ifftshift(ta2));

F=fftshift(fft(tu1.*tu2));

s=(4*m-1)*F(m+1:3*m-1);

end