function [abfull,abextended,abshort] = conv2cos(a,b)
% for cosine series in 2D
% abfull is just the full cosine convolution
% abextended is are all modes, not just the positive ones
% abshort has size(a) if size(a)==size(b)

Na = size(a)-1;
Nax=Na(1);
Nay=Na(2);
Nb = size(b)-1;
Nbx=Nb(1);
Nby=Nb(2);
Nx=max(Nax,Nbx);
Ny=max(Nay,Nby);
aa=altzeros([2*Nx+1,2*Ny+1],a(1));
bb=aa;

Sx=Nx+1;
Sy=Ny+1;

aa(Sx:Sx+Nax,Sy:Sy+Nay)=a;
aa(Sx-Nax:Sx,Sy:Sy+Nay)=flip(a,1);
aa(Sx:Sx+Nax,Sy-Nay:Sy)=flip(a,2);
aa(Sx-Nax:Sx,Sy-Nay:Sy)=flip(flip(a,1),2);

bb(Sx:Sx+Nbx,Sy:Sy+Nby)=b;
bb(Sx-Nbx:Sx,Sy:Sy+Nby)=flip(b,1);
bb(Sx:Sx+Nbx,Sy-Nby:Sy)=flip(b,2);
bb(Sx-Nbx:Sx,Sy-Nby:Sy)=flip(flip(b,1),2);


[~,abextended]=convtensor(aa,bb);

S2x=2*Nx+1;
S2y=2*Ny+1;
abextended=real(abextended);
abfull=abextended(S2x:S2x+Nax+Nb,S2y:S2y+Nay+Nb);
if Nax==Nbx && Nay==Nby
    abshort=abfull(1:Nax+1,1:Nay+1);
end

return
