function [ act, PLS, PX ] = nk_PLS_config(NM, PLS, PX, defaultsfl, parentstr)

algostr = 'pls';
UseLabel = true;
V = nk_MakeDummyVariables(NM.label);
cu = 1; cv = 1;
act=0;

if ~defaultsfl
    
    if isempty(PLS), [~, PLS ] = nk_PLS_config(NM, [], [], true); end
    algostr = PLS.algostr;
    cu = PLS.cu;
    cv = PLS.cv;
    UseLabel = PLS.uselabel;
    uselabels = {'NM labels','NM covariate matrix', 'User-defined matrix'};
    V = PLS.V;
    
    mn_str = sprintf('Choose PLS algorithm [ %s ]', algostr); mn_act = 1;
    mn_str = sprintf('%s|Define V matrix for covariance computation[ %s ]', mn_str, uselabels{UseLabel}); mn_act = [mn_act 2];
    if UseLabel == 2
        if isfield(NM,'covsel') && ~isempty(NM.covsel)
            colselstr = sprintf('%g indices selected',numel(NM.covsel));
        else
            colselstr = 'undefined';
        end
        mn_str = sprintf('%s|Provide column indices [ %s ]', mn_str, colselstr); mn_act = [mn_act 3];
    elseif UseLabel == 3
        if isfield(NM,'V') && ~isempty(NM.V)
            matselstr = sprintf('%gx%g matrix provided ',size(NM.V,1), size(NM.V,2));
        else
            matselstr = 'undefined';
        end
        mn_str = sprintf('%s|Provide user-defined matrix [ %s ]', mn_str, matselstr); mn_act = [mn_act 4];
    end
    
 	if strcmp(algostr,'spls')
        mn_str = sprintf('%s|Choose parameter(s) for sparsity constraint on U [ %s ]', mn_str, nk_ConcatParamstr(cu) ); mn_act = [mn_act 5];
        mn_str = sprintf('%s|Choose parameter(s) for sparsity constraint on V [ %s ]', mn_str, nk_ConcatParamstr(cv) ); mn_act = [mn_act 6];
    end

    nk_PrintLogo
    mestr = 'Partial least squares configuration'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', mn_str, mn_act);
    switch act
        case 1
            PLS.algostr = char(nk_input('Choose PLS method',0,'m','PLS|Sparse PLS',{'pls','spls'}));
        case 2
            if isfield(NM,'covars') && ~isempty(NM.covars)
                labeldesc = 'NM labels|NM covariate matrix|User-defined matrix'; labelnum = 1:3;
            else
                labeldesc = 'NM labels|User-defined matrix'; labelnum = 1:2;
            end
            PLS.uselabel = nk_input('Choose behavioral matrix for PLS',0,'m',labeldesc, labelnum, UseLabel);
            if PLS.uselabel == 1, PLS.V = nk_MakeDummyVariables(NM.label); end
        case 3
            if isfield(PLS,'covsel'), covsel = PLS.covsel; else, covsel = 1:size(NM.covars,2); end
            PLS.covsel = nk_SelectCovariateIndex(NM, covsel, 1);
            PLS.V = NM.covars(:,PLS.covsel);
        case 4
            PLS.V = nk_input('Provide a side label matrix for covariance matrix computation',0,'e',V,[size(NM.label,1) Inf]); 
        case 5
            PLS.cu = nk_input('Choose parameter (range) for sparsity constraint on U',0,'e',cu);
        case 6
            PLS.cv = nk_input('Choose parameter (range) for sparsity constraint on V',0,'e',cv);
           
    end
else
    PLS.algostr = algostr;
    PLS.uselabel = UseLabel;
    PLS.V = V;
    PLS.cu = cu;
    PLS.cv = cv;
    PX = nk_AddParam([], [], [], PX,'reset');
end

if strcmp(PLS.algostr,'spls')
    PX = nk_AddParam(PLS.cu,'SPLS-cu', 1, PX); 
    PX = nk_AddParam(PLS.cv,'SPLS-cv', 1, PX); 
else
    PX = nk_AddParam([],[],[],PX,'reset');
end