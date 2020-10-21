function [ThreshPerc, ThreshProb, Time] = EvalCoxPH(Models)

[ix, jx, nclass] = size(Models);
Time = cell(ix,jx,nclass);
ThreshPerc = zeros(ix,jx,nclass);
ThreshProb = zeros(ix,jx,nclass);
for h=1:nclass
    for k = 1:ix
        for l= 1:jx
            ThreshPerc(k,l,h) = Models{k,l,h}{1}.cutoff.val; 
            ThreshProb(k,l,h) = Models{k,l,h}{1}.cutoff.prob; 
            Time{k,l,h} = Models{k,l,h}{1}.predicted_time; 
        end
    end
end