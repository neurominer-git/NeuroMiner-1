function sY = nk_CopyY2Struct(Y, cv,i,j,kbin, sY, labels)

binmode=1;
TrInd = cv.TrainInd{i,j};
TsInd = cv.TestInd{i,j};

[iy,jy] = size(cv.cvin{i,j}.TrainInd);
try
    cpus = str2double(getenv('OMP_NUM_THREADS'));
catch
    cpus =1;
end

for u=1:kbin
    
   sY{i,j,u}.Tr = cell(iy,jy);
   sY{i,j,u}.CV = cell(iy,jy);
   sY{i,j,u}.Ts = cell(iy,jy);
   sY{i,j,u}.TrL = cell(iy,jy);
   sY{i,j,u}.CVL = cell(iy,jy);
   sY{i,j,u}.TsL = cell(iy,jy);
   
   for k=1:iy
       
       for l=1:jy
           
           switch binmode 
               
               case 0
                    fprintf('\n\t\tCV1 partition [perm = %g, fold = %g]', k,l)
                    TrInInd = TrInd(cv.cvin{i,j}.TrainInd{k,l});
                    CVInInd = TrInd(cv.cvin{i,j}.TestInd{k,l}); 

                    sY{i,j}.Tr = Y(TrInInd,:);  
                    sY{i,j}.CV = Y(CVInInd,:);
                    sY{i,j}.Ts = Y(TsInd,:);
                    sY{i,j}.TrL{k,l} = labels(TrInInd);
                    sY{i,j}.CVL{k,l} = labels(CVInInd);
                    sY{i,j}.TsL{k,l} = labels(TsInd);
                   
               case 1
                    fprintf('\n\t\tInner CV partition [perm = %g, fold = %g], bin. class.: %s.', ...
                       k,l, cv.class{i,j}{u}.groupdesc)
                    TrInBinInd = TrInd(cv.class{i,j}{u}.TrainInd{k,l});
                    CVInBinInd = TrInd(cv.class{i,j}{u}.TestInd{k,l});
                    TsInBinInd = TsInd(cv.classnew{i,j}{u}.ind);
                    
                    %sY{i,j,u}.Tr{k,l} = Y(TrInBinInd,:);  
                    %sY{i,j,u}.CV{k,l} = Y(CVInBinInd,:);
                    %sY{i,j,u}.Ts{k,l} = Y(TsInBinInd,:);
                    sY{i,j,u}.Tr{k,l} = extrX(Y,TrInBinInd,cpus);  
                    sY{i,j,u}.CV{k,l} = extrX(Y,CVInBinInd,cpus);
                    sY{i,j,u}.Ts{k,l} = extrX(Y,TsInBinInd,cpus);
                    
                    sY{i,j,u}.TrL{k,l} = cv.class{i,j}{u}.TrainLabel{k,l};
                    sY{i,j,u}.CVL{k,l} = cv.class{i,j}{u}.TestLabel{k,l};
                    sY{i,j,u}.TsL{k,l} = cv.classnew{i,j}{u}.label;
                   
           end
       end
   end
end

return