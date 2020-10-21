function OUT  = nk_SetFusionMode2(dat, analysis, F, nF, curmodal, oocvind)
global FUSION PREPROC VIS STACKING

if ~exist('oocvind','var'), oocvind = []; end; varstr = [];
if ~exist('curmodal','var') || isempty(curmodal), curmodal = 1; end
stk_flag = false; if ~isempty(STACKING) && STACKING.flag == 1, stk_flag = true; end
OUT.curmodal = curmodal;

%% Get Training / CV data (Y) & Build modality suffix
switch FUSION.flag
    
    case {0,3}
    
        varstr = ['_var' num2str(F(curmodal))]; tF = F(curmodal);
        disp_str = sprintf('\nPROCESSING OF MODALITY #%g', tF);
        if isfield(analysis,'GDdims'),
            switch FUSION.flag
                case 0
                     OUT.analysis = analysis.GDdims{1}; 
                case 3
                    if tF>numel(analysis.GDdims)
                        fprintf('\n'); warning('Requested analysis of modality %g which does not exist',tF)
                        OUT.analysis = [];
                    else
                        OUT.analysis = analysis.GDdims{tF}; 
                    end
            end
        end

    case 1
        
        disp_str = sprintf('\nPROCESSING OF CONCATENATED MODALITIES: EARLY FUSION -> '); 
        if numel(F)>3
            varstr = sprintf('_var_concat-%g',numel(F));
        else
            for i=1:numel(F)
                varstr = [varstr '_var' num2str(F(i))]; fprintf('#%g', F(i)); 
            end
        end
        tF = F;
        if isfield(analysis,'GDdims'), OUT.analysis = analysis.GDdims{1}; end

    case 2
        if nF > 1
            tF = F(curmodal); varstr = ['_var' num2str(tF)]; 
            disp_str = sprintf('\nPROCESSING OF MODALITY #%g', tF);
        else
            tF = F;
            disp_str = sprintf('\nPROCESSING OF CONCATENATED MODALITIES: INTERMEDIATE FUSION ->');
            for i=1:numel(F)
                varstr = [varstr '_var' num2str(F(i))]; fprintf('#%g', F(i)); 
            end
        end
        if isfield(analysis,'GDdims'), OUT.analysis = analysis.GDdims{1}; end

end

% Check whether stacking has to be activated
if stk_flag, 
    OUT.analyses = dat.analysis(STACKING.sel_anal); 
    OUT.nD = 0;
    for i=1:numel(OUT.analyses)
        OUT.nD = OUT.nD + numel(OUT.analyses{i}.GDdims);
    end
    OUT.stacking = true;
    varstr = '_varM'; 
    fprintf('\nPROCESSING OF PREDICTIONS FOR STACKING');
else
    OUT.analyses = [];
    OUT.nD = [];
    OUT.stacking = false;
    fprintf('%s',disp_str);
end

% Get back data and data descriptors according to data fusion mode
[OUT.X, OUT.labels] = nk_CompatY2(dat, tF, oocvind, FUSION);

% Initialize globals based on data fusion mode
switch FUSION.flag
    case {0,1,3}
        % Here only the first modality params will be initialized
        nk_SetupGlobVars2(analysis.params, 'setup_strat', 0, tF(1)); 
    case 2
        nk_SetupGlobVars2(analysis.params, 'setup_strat', 0, tF); 
end

% Check whether stacking is active and allow user to compute 
% fixed feature-order relevance metrics
OUT.fixedOrder = nk_CheckFixedFeatOrderStackInput(OUT.stacking, OUT.analyses);

OUT.PREPROC     = PREPROC;
OUT.VIS         = VIS;
OUT.F           = F;
OUT.tF          = tF;
OUT.l           = length(OUT.labels);    % # of subjects
if isfield(dat,'covars'),
    OUT.covars      = dat.covars;            % Covariates    
else
    OUT.covars = [];
end
OUT.covars_oocv = [];
OUT.desc_oocv   = [];
if ~isempty(oocvind), 
    if isfield(dat.OOCV{oocvind},'covars')
        OUT.covars_oocv = dat.OOCV{oocvind}.covars;
    end
    OUT.desc_oocv = dat.OOCV{oocvind}.desc;
end
OUT.featnames   = [];

if stk_flag
    if isfield(STACKING,'featnames'), 
        try 
            OUT.featnames{1} = STACKING.featnames(STACKING.sel_anal); 
        catch
            OUT.featnames{1} = STACKING.featnames; 
        end
    end
else
    if isfield(dat,'featnames') && ~isempty(dat.featnames(OUT.tF)) 
        OUT.featnames = dat.featnames(OUT.tF); 
    end
end

OUT.varstr = varstr;

