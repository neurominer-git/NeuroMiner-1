function Out = nk_CompSitePerf(G,L,P,func)

Out = zeros(size(G,2),1);
for i=1:size(G,2)
    Li = L(G(:,i)==1);
    Pi = P(G(:,i)==1);
    Out(i) = feval(func,Li,Pi);
end