function cout = ChebDerCoeffs(c,rigorous)
    [n, m] = size(c);
    if rigorous>0
      cout = intval(zeros(n-1, m));
      w = intval(repmat(2*(1:n-1)', 1, m));
    else
      cout = (zeros(n-1, m));
      w = repmat(2*(1:n-1)', 1, m);
    end
    v = w.*c(2:end,:);
    cout(n-1:-2:1,:) = vcumsum(v(n-1:-2:1,:));
    cout(n-2:-2:1,:) = vcumsum(v(n-2:-2:1,:));
    cout(1,:) = .5*cout(1,:);
end

function csum = vcumsum(V)
n=size(V,1);
csum=tril(ones(n))*V;
end

