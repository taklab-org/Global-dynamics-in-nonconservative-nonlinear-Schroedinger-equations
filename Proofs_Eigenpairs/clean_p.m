function p = clean_p(p)

S = (abs(p)<=1e-15);

p(S)=0;

end

