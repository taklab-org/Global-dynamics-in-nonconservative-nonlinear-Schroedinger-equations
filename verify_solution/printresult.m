sizeofy = size(y,1);
for i=1:sizeofy
  fprintf("%d & %d & %d & %.2f & %.2f & %.2f & %.2g & %.5g & %.5g & %.5g\\\\\n",...
    y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),y(i,6),y(i,7),y(i,8),y(i,9),y(i,10))
end