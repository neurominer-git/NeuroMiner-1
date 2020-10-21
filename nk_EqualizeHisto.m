function [tI, I] = nk_EqualizeHisto(Eq, V, I, modeflag)


if ~exist('V','var') || isempty(V), error('target label vector missing or undefined'); end

switch modeflag
    
    case 'regression'

        sV = sort(V,'ascend');

        for i=1:numel(sV)
           X = numel(find( V > sV(i)));
           if X <= Eq.MinCount, break, end
        end
        OptMaxV = sV(i-1);

        for i=numel(sV):-1:2
           X = numel(find( V < sV(i)));
           if X <= Eq.MaxCount, break, end
        end
        OptMinV = sV(i+1);

        Edges = OptMinV : (OptMaxV - OptMinV) / Eq.BinCount : OptMaxV;
        nE = numel(Edges);
        N = histc( V, Edges) ; N(end) = [];
        % fprintf('\nHistogram analysis of target labels:')
        % for i=1:numel(N)-1
        %     fprintf('\n%1.1f<->%1.1f =>\t%g', Edges(i), Edges(i+1), N(i));
        % end
        minN = min(N); tI = [];
        for i = 1:nE-1
            if i == nE
                ind = find(V >= Edges(i) & V <= Edges(i+1));
            else
                ind = find(V >= Edges(i) & V < Edges(i+1));
            end
           if N(i) > minN
               rInd = randperm(numel(ind));
               rInd = rInd(1:minN);
               tI = [ tI; I(ind(rInd)) ] ;
           else
               tI = [ tI; I(ind) ] ;
           end
        end
        ind = randperm(numel(tI));
        tI = tI(ind);
        I = setdiff(I, tI);
    case 'classification'
        
        uL = unique(V); nuL = numel(uL); suL= zeros(1,nuL); ruL = suL;
        for i = 1:nuL, suL(i) = sum(V==uL(i)); end
        [~,indmin] = min(suL); rI = [];
        targsize = floor(suL(indmin)*Eq.posnegrat);
        for i = 1:nuL, 
            ruL(i) = suL(i)/suL(indmin);
            if i==indmin || ruL(i)<Eq.posnegrat, continue; end
            remsize = floor(suL(i) - targsize);
            f = find(V==uL(i)); remind = f(randperm(numel(f),remsize));
            rI = [rI; remind];
        end
        tI = I(rI);
        I(rI) = [];
end
