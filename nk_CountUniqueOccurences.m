function [ID, NUM, TBL] = nk_CountUniqueOccurences(Y,PTH)
% =========================================================================
% function [ID, NUM, TBL] = nk_CountUniqueOccurences(Y,PTH)
% -------------------------------------------------------------------------
% Count unique occurences in Y and return labels as ID, along with their 
% frequencies (NUM). If PTH is provided write out frequency table to disk.
% =========================================================================
% (c) Nikolaos Koutsouleris, 01/2017
ID = unique(Y); nU=size(ID,1);
NUM = zeros(nU,1);
for i = 1:nU, NUM(i) = sum(Y==ID(i)); end
[~,ind] = sort(NUM,'descend');
ID = ID(ind); NUM=NUM(ind);
if exist('PTH','var')
   TBL = table(ID, NUM);
   P = fileparts(PTH);
   if exist(P,'dir')
        writetable(TBL,PTH)
   end
else
    TBL=[];
end