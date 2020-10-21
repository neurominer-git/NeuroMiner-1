function sparsity=computeSparsity(Y)
% compute sparsity

absY=abs(Y);
clear('Y');
maxY=max(max(absY));
[r,c]=size(absY);
sparsity=sum(sum(absY<=0.05*maxY))/(r*c);
end