function [CritGain, ExamFreq, PercThreshU, PercThreshL, AbsThreshU, AbsThreshL] = EvalSeqOpt(Models)

[ix, jx, nclass] = size(Models);
ExamFreq = cell(ix,jx,nclass);
CritGain = ExamFreq;
PercThreshU = ExamFreq;
PercThreshL = ExamFreq;
AbsThreshU = ExamFreq;
AbsThreshL = ExamFreq;

for h=1:nclass
    for k = 1:ix
        for l= 1:jx
            CritGain{k,l,h} = Models{k,l,h}{1}.allOPTs; 
            ExamFreq{k,l,h} = Models{k,l,h}{1}.examsfreq; 
            PercThreshU{k,l,h} = Models{k,l,h}{1}.optuvec; 
            PercThreshL{k,l,h} = Models{k,l,h}{1}.optlvec;
            AbsThreshU{k,l,h} = Models{k,l,h}{1}.optuthr;
            AbsThreshL{k,l,h} = Models{k,l,h}{1}.optlthr;
        end
    end
end