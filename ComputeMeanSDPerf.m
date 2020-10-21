function R = ComputeMeanSDPerf(NM, EXT, modind, PERFCRIT)

R.PERFs=[]; cnt=1;
if ~isempty(EXT)
    L = EXT.L;
    P = EXT.P;
    ind = EXT.ind;
    nInd = numel(unique(ind));
else
    L = NM.label; L(L==2)=-1;
    ind = NM.TrainParam.RAND.CV2LCO.ind;
    nInd = numel(unique(ind));
    P = NM.analysis{modind}.GDdims{1}.BinClass{1}.mean_predictions;
end
switch PERFCRIT
        case {1, 3, 4, 5, 'SENSITIVITY', 'BAC', 'PSI','AUC'}
            COMPARATOR = 1;
            switch PERFCRIT
                case 1
                    PERFCRIT='SENSITIVITY';
                case 3
                    PERFCRIT='BAC';
                case 4
                    PERFCRIT='AUC';
                case 5
                    PERFCRIT='PSI';
            end
        case {2,'SPECIFICITY'}
            COMPARATOR = -1;
            switch PERFCRIT
                case 2
                    PERFCRIT='SPECIFICITY';
            end
end
for i=1:nInd
    indi = ind==i;
    if any(L(indi)==COMPARATOR), 
        R.PERFs(cnt) = feval(PERFCRIT,L(indi), P(indi));
        cnt=cnt+1;
    end 
end
R.mean      = mean(R.PERFs);
R.std       = std(R.PERFs);
R.method    = PERFCRIT;
R.nInd      = nInd;
R.nSkipped  = nInd - numel(R.PERFs);
