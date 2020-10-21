%==========================================================================
%FORMAT OUT = MultiDichoFoldPerm(IN, OUT, F)
%==========================================================================
%MULTI-GROUP CLASSIFICATION PERFORMANCE EVALUATION                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(c) Nikolaos Koutsouleris, 07/2011

function OUT = MultiDichoFoldPerm(IN, OUT, F, LoopParam)
global VERBOSE MULTILABEL


if ~isfield(OUT,'mtr'), OUT.mtr         = cell(IN.nperms, IN.nfolds); end
if ~isfield(OUT,'mts'), OUT.mts         = cell(IN.nperms, IN.nfolds); end
if ~isfield(OUT,'mCVPred'), OUT.mCVPred = cell(IN.nperms, IN.nfolds); end
if ~isfield(OUT,'mTrPred'), OUT.mTrPred = cell(IN.nperms, IN.nfolds); end

nDicho          = zeros(IN.nclass,1); 
nTrSubj         = zeros(IN.nclass,1);
nCVSubj         = zeros(IN.nclass,1);
strout          = 'Multi-group classification performance';
Tr              = cell(IN.nclass,1);
CV              = cell(IN.nclass,1);
ngroups         = IN.ngroups; 

if exist('LoopParam','var') && ~isempty(LoopParam)
    PermVec = LoopParam.PermVec;
    FoldVec = LoopParam.FoldVec;
else
    PermVec = 1:IN.nperms;
    FoldVec = 1:IN.nfolds; 
end

PermNum = numel(PermVec);
FoldNum = numel(FoldVec);

for ii=1:PermNum % Loop through CV1 permutations

    for jj=1:FoldNum % Loop through CV1 folds
        
        i = PermVec(ii); j= FoldVec(jj);
        
        for curclass=1:IN.nclass % Loop through dichotomizers          
            % Optionally, select only some predefined subspaces
            if exist('F','var') && ~isempty(F)
                Fx = logical(F{i,j,curclass}); 
            else
                Fx = true(size(OUT.Trdecs{i,j,curclass},2),1); 
            end
            Tr{curclass} = OUT.Trdecs{i,j,curclass}(:,Fx);
            CV{curclass} = OUT.CVdecs{i,j,curclass}(:,Fx);
            [nTrSubj(curclass) nDicho(curclass)] = size(Tr{curclass});
            nCVSubj(curclass) = size(CV{curclass},1);
        end
        
        mxNFea = max(nDicho); 
        OUT.mtr{i,j} = zeros(mxNFea,1); 
        OUT.mts{i,j} = zeros(mxNFea,1);
        OUT.mCVPred{i,j} = zeros(nCVSubj(1), mxNFea);
        OUT.mTrPred{i,j} = zeros(nTrSubj(1), mxNFea);
        
        CVL = IN.Y.mCVL{i,j}(:,MULTILABEL.curdim); TrL = IN.Y.mTrL{i,j}(:,MULTILABEL.curdim);
        CVdecsCat = zeros(nCVSubj(1), IN.nclass);
        TrdecsCat = zeros(nTrSubj(1), IN.nclass);
        Classes = [];
        for k=1:mxNFea % Loop through feature subspaces
            
            indMx = k <= nDicho;
            
            % Construct multi-group decision matrix
            for curclass=1:IN.nclass
                if ~indMx(curclass), indFea = nDicho(curclass); else indFea = k; end
                TrdecsCat(:,curclass)   = Tr{curclass}(:,indFea);
                CVdecsCat(:,curclass)   = CV{curclass}(:,indFea);
                Classes = [Classes repmat(curclass,numel(indFea))];
            end
            
            % Evaluate multi-group performace for current feature subspace 
            % using specified method
            [OUT.mts{i,j}(k), OUT.mCVPred{i,j}(:,k)] = nk_MultiEnsPerf(CVdecsCat, sign(CVdecsCat), CVL, Classes, ngroups);
            [OUT.mtr{i,j}(k), OUT.mTrPred{i,j}(:,k)] = nk_MultiEnsPerf(TrdecsCat, sign(TrdecsCat), TrL, Classes, ngroups);
            
            if VERBOSE, fprintf('\n%s => CV1 [%g, %g, %g/%g subspaces]: %1.2f', strout, i, j, k, mxNFea, OUT.mts{i,j}(k)), end
        end
    end
end

end