function norm_a = dual_norm(a,inu)

m = length(a);

nu_power = [1 2*inu.^abs(1:m-1)];

norm_a = intval(max(sup(a./nu_power)));

end

