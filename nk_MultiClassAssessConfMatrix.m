function performance = nk_MultiClassAssessConfMatrix(confmatrix, label, pred, errs, crit)

if ~exist('crit','var') || isempty(crit),
    crit='all';
end
sumconf=sum(sum(confmatrix));
P=0; inan = isnan(pred); pred(inan) = []; label(inan)=[];
% Calculate sub confusion matrix for EACH GROUP vs the REST
uL = unique(label);
for i=1:size(confmatrix,1)
    labeli = label; labeli(label~=uL(i)) = 0;labeli(label==uL(i)) = 1;
    predi = pred; predi(pred~=uL(i)) = 0; predi(pred==uL(i)) = 1;
	h = confmatrix(i,:);
	h(i) = []; FN = sum(h);
	TP = confmatrix(i,i);
	FP = sum(confmatrix(:,i)) - TP;
	TN = sumconf-TP-FP-FN;
	performance.class{i}.TP     = TP;
	performance.class{i}.TN     = TN;
	performance.class{i}.FP     = FP;
	performance.class{i}.FN     = FN;
    labeli(~labeli)=-1; predi(~predi)=-1;
    switch crit
        case 'all'
            performance.class{i} = ALLPARAM(labeli,predi);
        otherwise
            performance.class{i}.(crit) = feval(crit,labeli,predi);
    end
    P = P + TP;
end
performance.confusion_matrix    = confmatrix;
performance.errors              = errs;
performance.accuracy            = P*100/sumconf;
if any(strcmp({'all','BAC'},crit)),
    performance.BAC                 = [];
    for i = 1:numel(performance.class)
        performance.BAC = [ performance.BAC performance.class{i}.BAC];
    end
    performance.BAC_SD = nanstd(performance.BAC);
    performance.BAC_Mean = nanmean(performance.BAC);
end

return