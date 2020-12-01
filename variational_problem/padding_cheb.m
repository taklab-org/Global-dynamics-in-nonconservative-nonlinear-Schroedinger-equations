function c = padding_cheb(c,n,N,n_pad_cheb)

c = reshape(c,n,2*N+1); 
c = [c;zeros(n_pad_cheb,2*N+1)];
c = reshape(c,(n+n_pad_cheb)*(2*N+1),1);

end

