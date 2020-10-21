function [oECOC, tECOC] = nk_OneVsOne(Dichotomizers, number_classes)

nc = size(Dichotomizers,2);
mx = number_classes*(number_classes-1)/2;
% This is the observed matrix
oECOC=zeros([number_classes nc]);
% This is the template matrix
tECOC=zeros([number_classes mx]);
counter=1;

% Create multi-group coding matrix (template) for one-vs-one
for i=1:number_classes-1
    for j=i+1:number_classes
        tECOC(i,counter)=1;
        tECOC(j,counter)=-1;
        counter=counter+1;
    end
end

% Assign dichotomizers to class vector
for i=1:mx
    ind = Dichotomizers == i; l = sum(ind);
    oECOC(:,ind) = repmat(tECOC(:,i),1,l);
end

return
