function OOCVres = nk_OOCVMeta(OOCVres, inp)
global META MULTILABEL CV MODEFL

if isempty(META) || ( ~META.flag || inp.batchflag )
    if inp.nF>1, fprintf('\nSkipping bagging step.'); end
    return
end
%aggmode = 1;
%nG = inp.nF;
nclass = inp.nclass;
Labels = double.empty(0,nclass); MultiLabels = [];

if isfield(inp,'labelOOCV') && ~isempty(inp.labelOOCV)
    switch MODEFL
        case 'classification'
            lu = [1 -1];
            Labels = nan(numel(inp.labelOOCV),inp.nclass);
            for curclass = 1:nclass
                cvu = CV.class{1,1}{curclass};
                for y = 1:numel(cvu.groups)
                    indy = inp.labelOOCV == cvu.groups(y); Labels(indy,curclass) = lu(y);
                end
            end
            MultiLabels = inp.labelOOCV;
        case 'regression'
            Labels = inp.labelOOCV;
    end
end
        
for j=1:MULTILABEL.dim
    for curclass = 1:nclass
        % The OOCV predictions are already stored in a bagged ensemble
        % we just copy them over to the META structure array
        L = Labels(:,curclass);
        P = mat2cell(OOCVres.predictions{curclass,j},ones(inp.nOOCVsubj,1), size( OOCVres.predictions{curclass, j}, 2) );
        switch MODEFL
            case 'classification'
                OOCVres.META{j}.BinClass{curclass} = nk_ComputeEnsembleProbability(P, L);
                if isfield(OOCVres,'multi_predictions')
                    mP = OOCVres.multi_predictions{j};
                    nGroups = unique(inp.labels); nGroups(isnan(nGroups))=[]; nGroups = numel(nGroups);
                    OOCVres.META{j} = nk_MultiPerfComp(OOCVres.META{j}, mP, MultiLabels, nGroups);
                end
            case 'regression'
                OOCVres.META{j}.Regr = nk_ComputeEnsembleProbability(P, L);
        end
    end
end