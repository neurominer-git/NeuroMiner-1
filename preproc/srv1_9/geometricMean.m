function [gm]=geometricMean(v)
% for example
% v=[1,2,3,NaN];
% [gm]=geometricMean(v)
% Yifeng Li
% July 04, 2012

len=numel(v);
lenNaN=len;
p=1; % product
for i=1:len
   if isnan(v(i))
      lenNaN=lenNaN-1;
      continue;
   end
   p=p*v(i);
end
if lenNaN>0
   gm=nthroot(p,lenNaN);
else
    gm=NaN;
end
end