function miY = nk_MI(Y, dat, SplitFlag)

[ix,jx] = size(dat.cv.TrainInd);
miY{ix,jx} = [];
MItype = dat.TrainParam.MI.type;

for i=1:ix

    for j=1:jx
    
        fprintf('\nScaling outer CV folds [%g,%g].',i,j)
        
        switch SplitFlag

            % Dimred. is performed on ALL data for each outer CV loop
            case 0 
                   TrInd = 1:size(Y,1);
                   CVInd=[];
                   
            % DimRed. is performed on the outer training fold data 
            % and an out-of-sample mapping is performed on the outer CV dataV loop
            case 1
                   TrInd = dat.cv.TrainInd{i,j};
                   CVInd= dat.cv.TestInd{i,j};
                   
            % Dimred. is performed on the inner training data and 
            % an out-of-sample mapping is performed on the inner and outer CV data. 
            case 2  
                   TrInd = dat.cv.TrainInd{i,j};
                   CVInd = dat.cv.TestInd{i,j};
                   [iy,jy] = size(dat.cv.cvin{i,j}.TrainInd);
                   miY{i,j}.inY{iy,jy}=[];
                   for k=1:iy
                       for l=1:jy
                           fprintf('\n\tScaling inner CV folds [%g,%g].',k,l)
                           TrInInd = dat.cv.TrainInd{i,j}(dat.cv.cvin{i,j}.TrainInd{k,l});
                           CVInInd = dat.cv.TrainInd{i,j}(dat.cv.cvin{i,j}.TestInd{k,l});
                           if iscell(Y)
                               tY = Y{i,j}.inY{k,l};
                           else
                               tY = Y;
                           end
                                switch MItype
                                    case 1
                                        miY{i,j}.inY{k,l} = nk_MIBag(tY(TrInInd), dat.label(TrInInd),0);
                                    case 2
                                        miY{i,j}.inY{k,l} = nk_MIPerm(tY(TrInInd), dat.label(TrInInd), 0, dat.TrainParam.MI.nperms);
                                end
                           end
                   end
        end
        if iscell(Y)
            tY=Y{i,j}.Y;
        else
            tY=Y;
        end
         switch MItype
            case 1
                miY{i,j}.Y = nk_MIBag(tY(TrInd), dat.label(TrInd), 0);
            case 2
                miY{i,j}.Y = nk_MIPerm(tY(TrInd), dat.label(TrInd), dat.TrainParam.MI.nperms);
         end
    end
end

return