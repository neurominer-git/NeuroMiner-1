function res = crossanalyse(MODEFL, L, analyses)

switch MODEFL
    case 'classification'
        ResField = 'BinClass';
    case 'regression'
        ResField = 'Regr';
end

nA = numel(analyses);
nclass = numel(analyses{1}.GDdims{1}.(ResField));
nsubj = size(L,1);
res.D = cell(1,nclass);
for curclass = 1:nclass
    res.D{curclass} = zeros(nsubj,nA);
    for a = 1:nA
        nG = numel(analyses{a}.GDdims); tD = zeros(nsubj, nG);
        if nG>1
            for b = 1:nG
                tD(:,b) = analyses{a}.GDdims{b}.(ResField){curclass}.mean_predictions; 
            end
            tD = mean(tD,[],2);
        else
            tD = analyses{a}.GDdims{1}.(ResField){curclass}.mean_predictions; 
        end
        res.D{curclass}(:,a) = tD;
   end
   cnt = 1;
   nI = sum(isnan(res.D{curclass}),2)>0; 
   res.D{curclass}(nI,:)=[];
   
   for ai=1:nA-1
       for aj = ai+1:nA
            res.anal(cnt).comparison = sprintf('Analysis %s vs. Analysis %s', analyses{ai}.id,analyses{aj}.id);
            [res.anal(cnt).pearson_r, res.anal(cnt).pearson_p ] = corrcoef(res.D{curclass}(:,ai), res.D{curclass}(:,aj));
            res.anal(cnt).mcnem2x2.mat =[ sum(res.D{curclass}(:,ai)>0 & res.D{curclass}(:,aj)>0) sum(res.D{curclass}(:,ai)>0 & res.D{curclass}(:,aj)<0); ...
                            sum(res.D{curclass}(:,ai)<0 & res.D{curclass}(:,aj)>0) sum(res.D{curclass}(:,ai)<0 & res.D{curclass}(:,aj)<0) ];  
            [res.anal(cnt).mcnem2x2.chi2crit, ...
                res.anal(cnt).mcnem2x2.chi2val, ... 
                res.anal(cnt).mcnem2x2.pval, ...
                res.anal(cnt).mcnem2x2.power ] = mcnemar( res.anal(cnt).mcnem2x2.mat);        
            res.anal(cnt).kappa = KAPPA(res.anal(cnt).mcnem2x2.mat(1,1), res.anal(cnt).mcnem2x2.mat(2,2), res.anal(cnt).mcnem2x2.mat(1,2), res.anal(cnt).mcnem2x2.mat(2,1));
            cnt=cnt+1;
       end
   end
    
end



