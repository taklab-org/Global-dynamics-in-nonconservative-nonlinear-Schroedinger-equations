sizeofy = size(y,1);
for i=1:sizeofy
  fprintf("%.2g & %.2g & %.5g & %.5g & %.5g & %.2f & %.2f & %.5g & %.5g & %.5g\\\\\n",...
    mid(y(i,1)),mid(y(i,2)),mid(y(i,3)),sup(y(i,4)),sup(y(i,5)),sup(y(i,6)),sup(y(i,7)),sup(y(i,8)),sup(y(i,9)),sup(y(i,10)))
end