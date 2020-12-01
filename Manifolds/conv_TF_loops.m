function ab = conv_TF_loops(a,b)

Size = size(a)-1; N = Size(1); M = Size(2);

ab = altzeros(size(a),a(1,1));

for k=0:N
    for m=0:M
        s_km = 0;
        for ell=0:m
            for k1=-N:N
                if abs(k-k1)<=N
                    s_km = s_km + a(abs(k1)+1,ell+1)*b(abs(k-k1)+1,m-ell+1);
                end
            end
        end
        
        ab(k+1,m+1) = s_km;
    end
end

end