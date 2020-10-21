function oECOC = nk_OneVsAll(Classes, number_classes)

nc = size(Classes,2);
oECOC=zeros([number_classes nc]);
tECOC=-1*ones([number_classes number_classes]);

% Create multi-group coding matrix for one-vs-one
for i=1:number_classes
    tECOC(i,i)=1;
end

% Assign dichotomizers to class vector
for i=1:number_classes
    ind = Classes == i; l = sum(ind);
    oECOC(:,ind) = repmat(tECOC(:,i),1,l);
end

return
