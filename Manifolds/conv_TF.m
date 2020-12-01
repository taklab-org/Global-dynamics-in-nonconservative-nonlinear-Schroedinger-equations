function [abfull,abshort] = conv_TF(a,b)
% for Taylor-Fourier series
% abfull is just the full convolution Taylor-Fourier product
% abshort has size(a) if size(a)==size(b)

Na = size(a)-1;
Nax = Na(1);
Nay = Na(2);
Nb = size(b)-1;
Nbx = Nb(1);
Nby = Nb(2);
Nx = max(Nax,Nbx);
Ny = max(Nay,Nby);
aa = altzeros([2*Nx+1,2*Ny+1],a(1));
bb = aa;

Sx = Nx+1;
Sy = Ny+1;

aa(1:Sx+Nax,Sy:Sy+Nay) = [flip(a(2:end,:),1);a]; 
bb(1:Sx+Nbx,Sy:Sy+Nby) = [flip(b(2:end,:),1);b];

[~,abextended] = convtensor(aa,bb);

S2x = 2*Nx+1;
S2y = 2*Ny+1;
%abextended = real(abextended);
abfull = abextended(S2x:S2x+Nax+Nbx,S2y:S2y+Nay+Nby);
if Nax==Nbx && Nay==Nby
    abshort=abfull(1:Nax+1,1:Nay+1);
end

return
