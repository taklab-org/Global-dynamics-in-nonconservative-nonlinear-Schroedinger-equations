function cout = ChebDerCoeffs_fourier(c)
    [n, m] = size(c);
    cout = (zeros(n-1, m));
    w = repmat(2*(1:n-1)', 1, m);
    v = w.*c(2:end,:);
    cout(n-1:-2:1,:) = vcumsum(v(n-1:-2:1,:));
    cout(n-2:-2:1,:) = vcumsum(v(n-2:-2:1,:));
    cout = [cout;zeros(1,m)];
end

function csum = vcumsum(V)
n=size(V,1);
csum=tril(ones(n))*V;
end

