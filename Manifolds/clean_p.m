function p = clean_p(p)

S = (abs(p)<=1e-16);

p(S)=0;

end

