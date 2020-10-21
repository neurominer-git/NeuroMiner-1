function [sY, params] = nk_Scale(Y, cv, i, j, kbin, params, labels)

binmode=1;
TrInd = cv.TrainInd{i,j};
TsInd = cv.TestInd{i,j};

if iscell(Y)
    sY=Y;
end

[iy,jy] = size(cv.cvin{i,j}.TrainInd);

for u=1:kbin
    
   sY{i,j,u}.Tr = cell(iy,jy);
   sY{i,j,u}.CV = cell(iy,jy);
   sY{i,j,u}.Ts = cell(iy,jy);
   params{i,j,u}.minTr = cell(iy,jy);
   params{i,j,u}.maxTr = cell(iy,jy);
   
   if ~isfield(Y{i,j,u}, 'TrL')
       sY{i,j,u}.TrL = cell(iy,jy);
       sY{i,j,u}.CVL = cell(iy,jy);
       sY{i,j,u}.TsL = cell(iy,jy);
   end
   
   for k=1:iy
       
       for l=1:jy
           
           switch binmode 
               
               case 0
                   
                   TrInInd = TrInd(cv.cvin{i,j}.TrainInd{k,l});
                   CVInInd = TrInd(cv.cvin{i,j}.TestInd{k,l}); 
                   if iscell(Y)
                       Tr = Y{i,j}.Tr{k,l};
                       CV = Y{i,j}.CV{k,l};
                       Ts = Y{i,j}.Ts{k,l};
                   else
                       Tr = Y(TrInInd,:);  
                       CV = Y(CVInInd,:);
                       Ts = Y(TsInd,:);
                   end
                   fprintf('\n\t\tCV1 partition [perm = %g, fold = %g]', ...
                           k,l)
                   [sY{i,j}.Tr{k,l}, sY{i,j}.CV{k,l}, sY{i,j}.Ts{k,l}, ...
                       params{i,j}.minTr{k,l}, params{i,j}.maxTr{k,l}] = nk_PerfScale(Tr, CV, Ts);

                   if ~isfield(sY{i,j},'TrL')
                        sY{i,j}.TrL{k,l} = labels(TrInInd);
                        sY{i,j}.CVL{k,l} = labels(CVInInd);
                        sY{i,j}.TsL{k,l} = labels(TsInd);
                   end
                   
               case 1

                   TrInBinInd = TrInd(cv.class{i,j}{u}.TrainInd{k,l});
                   CVInBinInd = TrInd(cv.class{i,j}{u}.TestInd{k,l});
                   TsInBinInd = TsInd(cv.classnew{i,j}{u}.ind);
                   
                   if iscell(Y)
                       Tr = Y{i,j,u}.Tr{k,l};
                       CV = Y{i,j,u}.CV{k,l};
                       Ts = Y{i,j,u}.Ts{k,l};                                         
                   else
                       Tr = extrX(sY,TrInBinInd);  
                       CV = extrX(sY,CVInBinInd);
                       Ts = extrX(sY,TsInBinInd);
                   end
                   
                   fprintf('\n\t\tCV1 partition [perm = %g, fold = %g], bin. class.: %s.', ...
                       k,l, cv.class{i,j}{u}.groupdesc)

                   [sY{i,j,u}.Tr{k,l}, sY{i,j,u}.CV{k,l}, sY{i,j,u}.Ts{k,l}, ...
                       params{i,j,u}.minTr{k,l}, params{i,j,u}.maxTr{k,l}] = nk_PerfScale(Tr, CV, Ts);

                   if isempty(sY{i,j,u}.TrL{k,l})
                        sY{i,j,u}.TrL{k,l} = cv.class{i,j}{u}.TrainLabel{k,l};
                        sY{i,j,u}.CVL{k,l} = cv.class{i,j}{u}.TestLabel{k,l};
                        sY{i,j,u}.TsL{k,l} = cv.classnew{i,j}{u}.label;
                   end
           end
       end
   end
end

return