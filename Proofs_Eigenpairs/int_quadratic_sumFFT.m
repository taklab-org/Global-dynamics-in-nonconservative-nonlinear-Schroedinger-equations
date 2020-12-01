function [a1a2]=int_quadratic_sumFFT(a1,a2)
% the convolution term using FFT

M=length(a1);
a1=[flip(a1(2:end));a1];
a2=[flip(a2(2:end));a2];

i2=intval('2');
i4=intval('4');

% We make sure that the inputs are powers of 2 %
M1=(2^ceil(log(4*M-1)/log(2))-4*M)/2;   % Hence 4*M+2*M1 is a power of 2 %

ta1=[intval(zeros(M+M1+1,1));a1;intval(zeros(M+M1,1))];
ta2=[intval(zeros(M+M1+1,1));a2;intval(zeros(M+M1,1))];

tu1=verifyfft(ifftshift(ta1),-1);
tu2=verifyfft(ifftshift(ta2),-1);
F=fftshift(verifyfft(tu1.*tu2,1));

a1a2 = (i4*intval(M)+i2*intval(M1))*F(M1+2*M+1:M1+3*M,1);

end