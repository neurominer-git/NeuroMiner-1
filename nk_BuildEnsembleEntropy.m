function [ind, origperf, finalperf] = nk_BuildEnsembleEntropy(decs, targs, TsInd, TsIdx, TsLabels)

global SVM RFE

% Build prediction value matrix for entire ensemble
ct = length(targs);
[D,T,L] = nk_Vals2Ind(decs,targs,TsInd,TsIdx,TsLabels);

S=1:ct;
lbind = any(T(:,S),2);
hx = sign(sum(D(lbind,S),2));
lb = sign(sum(L(lbind,S),2));
origperf = get_test_param(SVM.GridParam,hx,lb);

ambiguity(1) = nk_Ambiguity(T(:,S)); 
performance(1) = origperf;
MaxPerf = origperf;
MaxAmbi = ambiguity(1);
fprintf(1,'\nOriginal ensemble size: %g',ct)
ind = 1:ct;flg=0;

% Maximize CV performance by backward elimination of hypotheses
% As soon as ambiguity starts dropping stop the process
while ~isempty(S)
 
    l=length(S);
    A = zeros(l,1);
    AC = A;
    
    while l > 0
        kS = S; kS(l)=[];
        A(l) = nk_Ambiguity(T(:,kS));
        lbind = any(T(:,kS),2);
        hx = sign(sum(T(lbind,kS),2));
        
        % Check if there are 50% decisions
        % and through a coin
        if any(hx==0) > 0
           hx0 = hx==0; 
           shx0 = sum(hx0);
           i0 = rand(shx0,1);
           i0(i0>0.5) = 1; i0(i0~=1)=-1;
           hx(hx0) = i0;
        end
        lb = sign(sum(L(lbind,kS),2));
        AC(l) = get_test_param(SVM.GridParam,hx,lb);
        l=l-1;
    end


	[param, elind] = max(A);
    
    	switch RFE.EnsMethod.Crit
		case 1 % Maximize ambiguity (entropy)
			if param > MaxAmbi, flg = 1; end;
		case 2 % Maximize performance
			if AC(elind) > MaxPerf, flg = 1; end
		case 3 % Maximize both
			if param > MaxAmbi && AC(elind) > MaxPerf, flg = 1; end
	end	
	
	if flg
		ind = S;
		ambiguity(end+1) = A(elind);
		performance(end+1) = AC(elind);
		MaxPerf = AC(elind);
		MaxAmbi = A(elind);
		flg=0;
	end
	S(elind)=[];
end

if RFE.disp
    figure(8)
    subplot(2,1,1); plot(1:length(ambiguity),ambiguity)
    title('Ambiguity')
    subplot(2,1,2); plot(1:length(performance),performance)
    title('Performance')
    xlabel('Ensemble classifier size')
end

finalperf = MaxPerf;
fprintf(1,'\nFinal ensemble size: %g',length(ind))
fprintf(1,'\nDiversity: %1.3f',ambiguity(end))
fprintf(1,'\nOrig/Final Performance: %1.1f/%1.1f',origperf,finalperf)

return