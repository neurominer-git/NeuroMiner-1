function [mPerf, sdPerf, Params, bars] = extract_subparam_performance(P, Perf, Pind, curclass, plotflag, ha)
nclass=1;
if iscell(P) 
    if ~exist('curclass','var') || isempty(curclass)
        nclass = numel(P); mPerf = cell(nclass,1); sdPerf = cell(nclass,1); Params = cell(nclass,1);
        for curclass = 1:nclass
            [mPerf{curclass}, sdPerf{curclass}, Params{curclass}] = extract_perf(P{curclass}(:,Pind), Perf(:,curclass));
        end
    else
        if size(Perf,2)==numel(P),
            [mPerf, sdPerf, Params] = extract_perf(P{curclass}(:,Pind), Perf(:,curclass));
        else
            [mPerf, sdPerf, Params] = extract_perf(P{curclass}(:,Pind), Perf);
        end
    end
else
    if size(Perf,2)>1
        nPerf = size(Perf,2);
        mPerf = cell(nPerf,1); sdPerf = cell(nPerf,1); Params = cell(nPerf,1);
        for i=1:nPerf
            [mPerf{i}, sdPerf{i}, Params{i}] = extract_perf(P(:,Pind), Perf(:,i));
        end
    else
        [mPerf, sdPerf, Params] = extract_perf(P(:,Pind), Perf);
    end
end

if exist('plotflag','var') || ~isempty(plotflag)
   if ~exist('ha','var') || isempty(ha), 
       figure('Name','NM Performance Analysis Plot: Parameter subspaces'); hold on
       ha = gca; 
   end
   if (iscell(P) && nclass > 1) || iscell(mPerf)
      mPx = cell2mat(mPerf)'; sPx = cell2mat(sdPerf)'; Parx = cell2mat(Params');    
   
   else
      mPx = mPerf; sPx = sdPerf; Parx = Params;    
   end
   cla(ha);
   bars = bar(ha, mPx); %errorbar(repmat( (1:size(Parx,1))',1,nclass),mPx,sPx,'LineStyle','none');
   ha.XLim = [0.5 size(Parx,1)+0.5]; ha.YLimMode = 'auto';
   ha.XTick = 1:size(Parx,1); ha.XTickLabel=num2str(Parx(:,1));

end

function [mPerf, sdPerf, uPi] = extract_perf(Pi, Perf)

uPi = unique(Pi);
nuPi = numel(uPi);
if nuPi==1, warning('Only one parameter found'); end
mPerf = arrayfun(@(i) mean(Perf(Pi == uPi(i))), 1:nuPi);
sdPerf = arrayfun(@(i) std(Perf(Pi == uPi(i))), 1:nuPi);