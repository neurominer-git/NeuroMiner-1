function [confmatrix, errs, pred] = nk_DecideClassMem(class,labels,indices,predictions)

global MULTI SVM

ml = max(labels); nclass = size(class,1); lx = length(labels); confmatrix = zeros(ml,ml);

switch MULTI.decisiontype

case 1 % One versus one maximum wins
    if strcmp(SVM.prog,'LIBSVM') && strcmp(SVM.Optimization.b,' -b 1')
        predictions(predictions == 0.5) = predictions(predictions == 0.5) + 0.001;
        predictions(predictions~=0) = predictions(predictions~=0) - 0.5;
    end
    [dum,maxI] = max(abs(predictions),[],2);
    ind0 = false(lx,nclass); ind0(sub2ind([lx nclass], 1:lx, maxI')) = true;
    indlabels=zeros(size(indices));
    
    sp = sign(predictions);
    for i=1:size(indices,2) % Loop through binary classifiers
        indlabels(sp(:,i) == 1,i) =  class{i}.groups(1);
        indlabels(sp(:,i) ==-1,i) =  class{i}.groups(2);
    end
    
    %errs = zeros(lx,1); pred = zeros(lx,1);
    errs = sum(indices.*ind0,2) ~= sum(sp.*ind0,2);
    pred = sum(indlabels.*ind0,2);
    
    for i=1:lx
%         errs(i) = indices(i,ind0(i,:)) ~= sp(i,ind0(i,:));
%         pred(i) = indlabels(i,ind0(i,:));
        confmatrix(labels(i),pred(i)) = confmatrix(labels(i),pred(i)) + 1;
    end
    
case 2 % Voting method
    
    predictionsbin = zeros(size(predictions));
    predictionsbin(predictions<0) = -1;
    predictionsbin(predictionsbin==0) = 1;
	
    % not completed
    
case 3 % DAG method
    
    classmem=zeros(size(labels));
    for i=1:length(labels)
        for j=1:length(class)
            if indices(i,j)==0, continue, end;
            if predictions(i,j) < 0
                classmem(i) = class{j}.groups(1);
            else
                classmem(i) = class{j}.groups(2);
            end
        end
        conf_x = labels(i);
        confmatrix(conf_x,classmem(i)) = confmatrix(conf_x,classmem(i)) + 1;
    end
    errs = classmem ~= labels;

end
