function OUT = DetrendFoldPerm(IN, OUT)

global SVM MODEFL

if isfield(SVM,'Post') && isfield(SVM.Post,'Detrend') && SVM.Post.Detrend
    OUT.detrend.flag = true;
else
    OUT.detrend.flag = false;
end

if OUT.detrend.flag

    PermVec = 1:IN.nperms;
    FoldVec = 1:IN.nfolds;
    ClassVec = 1:IN.nclass;

    PermNum = numel(PermVec);
    FoldNum = numel(FoldVec);
    ClassNum = numel(ClassVec);

    P = []; D = []; L = [];

    for ii=1:PermNum % Loop through CV1 permutations

        for jj=1:FoldNum % Loop through CV1 folds

            i = PermVec(ii); j= FoldVec(jj);

            for ccurclass=1:ClassNum % Loop through dichotomizers

                curclass = ClassVec(ccurclass);
                %CVInd = IN.Y.CVInd{i,j}{curclass};
                %TrInd = IN.Y.TrInd{i,j}{curclass};
                zu = size(OUT.Trtargs{i,j,curclass},2);

                for z = 1 : zu
                    
                    D = [ D; OUT.CVdecs{i,j,curclass}(:,z) ] ;
                    P = [ P; OUT.CVtargs{i,j,curclass}(:,z) ] ;
                    L = [ L; IN.Y.CVL{i,j}{curclass} ] ; 
                    
                end

            end
        end
    end
    
    switch MODEFL
        case 'regression'
            [~, OUT.detrend.beta, OUT.detrend.p] = nk_DetrendPredictions2([], [], P, L);
            
        case 'classification'
            [ ~,~, T, ~, ~, OUT.detrend.indexthresh] = perfcurve2(L,D,1);
            OUT.detrend.thresh = T(OUT.detrend.indexthresh);
            %Pcorr = P - OUT.detrend.thresh;
             for ii=1:PermNum % Loop through CV1 permutations
                for jj=1:FoldNum % Loop through CV1 folds
                    i = PermVec(ii); j= FoldVec(jj);
                    for ccurclass=1:ClassNum % Loop through dichotomizers
                        curclass = ClassVec(ccurclass);
                        zu = size(OUT.Trtargs{i,j,curclass},2);
                        for z = 1 : zu
                            OUT.Trdecs{i,j,curclass}(:,z) =  OUT.Trdecs{i,j,curclass}(:,z) - OUT.detrend.thresh ;
                            OUT.Trtargs{i,j,curclass}(:,z) = sign( OUT.Trdecs{i,j,curclass}(:,z) );
                            OUT.CVdecs{i,j,curclass}(:,z) =  OUT.CVdecs{i,j,curclass}(:,z) - OUT.detrend.thresh ;
                            OUT.CVtargs{i,j,curclass}(:,z) = sign( OUT.CVdecs{i,j,curclass}(:,z) );
                        end
                    end
                end
            end
    end

end

end