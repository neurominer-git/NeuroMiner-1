function TEMPL = nk_CreatePreprocTemplate(Y, label)
global CV RAND PREPROC MODEFL VERBOSE

if isfield(PREPROC,'TEMPLPROC') && ~isempty(PREPROC.TEMPLPROC) && PREPROC.TEMPLPROC
    if VERBOSE, fprintf('\nCreate full population preprocessing template for Procrustes aligmnent'); end
    nshelves = numel(Y);
    switch MODEFL
        case 'classification'
            SrcParam.MultiTrainLabel    = label;
        case 'regression'
            SrcParam.TrainLabel         = label;
    end
    if PREPROC.BINMOD
        ukbin = numel(CV.class{1,1}); SrcParam.binmult = 1;
    else
        ukbin = 1; SrcParam.binmult = 0;
    end
    SrcParam.CV1perm = 0; SrcParam.CV1fold = 0;
    TEMPL.Tr = cell(1,ukbin); TEMPL.Param = cell(1,ukbin);
    for curclass = 1 : ukbin
        SrcParam.u = curclass;
        ind1 = find(label == CV.class{1,1}{curclass}.groups(1)); ind2 = find(label == CV.class{1,1}{curclass}.groups(2));
        if RAND.Decompose ~=9
            SrcParam.BinaryTrainLabel = zeros(numel(ind1)+numel(ind2),1);
            SrcParam.BinaryTrainLabel(1:numel(ind1)) = 1; SrcParam.BinaryTrainLabel(numel(ind1)+1:end) = -1;
        end
        SrcParam.TrX = [ind1;ind2];
        for shelf = 1: nshelves
            InputParam.Tr{shelf} = Y{shelf}(SrcParam.TrX,:);
        end
        [TEMPL.Tr{curclass}, TEMPL.Param{curclass}] = nk_GenPreprocSequence(InputParam, PREPROC, SrcParam);
    end    
else
    TEMPL = [];fprintf('\nRun chain without template processing.')
end