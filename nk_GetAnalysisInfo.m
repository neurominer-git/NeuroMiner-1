function str = nk_GetAnalysisInfo(dat, analysis)

if iscell(analysis.params.datadescriptor), nvar = length(analysis.params.datadescriptor); else nvar = 1; end
showmodalmax = Inf;
showmodalvec = 1:nvar;
params = analysis.params;
FUSION = params.TrainParam.FUSION;
if isempty(FUSION), FUSION.flag = 0; end
nF=1; str=[];
if nvar > 1
    nF = numel(FUSION.M);
    if nF>1,
        strmod = sprintf('%g MODALITIES',nF);
    else
        strmod = 'MODALITY';
    end
    str{1} = sprintf('ANALYSIS OPERATES ON %s', strmod); 
    for j=1:nF, 
        if j > showmodalmax, fprintf(' ...'); break; end;
        str{end+1} = sprintf(' #%g', FUSION.M(j)); 
    end
    str{end+1} = fprintf('\n'); 
end
e = nk_GetParamDescription2(dat, params,'cv');

% Loop through variates
for jj=1:nF

    j = showmodalvec(jj);

    if FUSION.flag == 1 && j>1
        continue; 
    elseif FUSION.flag == 3
        params.TrainParam = analysis.params.TrainParam.STRAT{j};
    else
        params.TrainParam = analysis.params.TrainParam;
    end

    switch FUSION.flag
        case {0,2}
            % Get info about modality j
            d = nk_GetParamDescription2(dat, params, 'VarDesc', [], FUSION.M(j));
            d = nk_GetParamDescription2(dat, params.TrainParam.PREPROC{FUSION.M(j)},'PreProc', d, j);
        case 1
            % Get info about all modalities
            d = nk_GetParamDescription2(dat, params,'VarDesc', [], FUSION.M);
            d = nk_GetParamDescription2(dat, params.TrainParam.PREPROC{j},'PreProc', d, 1);
        case 3
            % Get info about modality j
            d = nk_GetParamDescription2(dat, params,'VarDesc', [], FUSION.M);
            d = nk_GetParamDescription2(dat, params.TrainParam.PREPROC,'PreProc', d);                    
    end

    if FUSION.flag == 1
        mxl = 0;
        for i=1:numel(d.datadescriptor)
            str{end+1} = sprintf('MODALITY %g : %s \n', FUSION.M(i), d.datadescriptor{i});
            mxli = size(str{end},2);
            if mxli > mxl, mxl = mxli; end
        end
    else
        if numel(d.datadescriptor) > 1,
            datdesc = d.datadescriptor{j};
        else
            datdesc = d.datadescriptor{1};
        end
        str{end+1} = sprintf('MODALITY %g : %s \n', FUSION.M(j), datdesc);
        mxl = size(str{end},2);
    end

    str{end+1} = sprintf('%s \n',repmat('*',1,mxl));
    str{end+1} = sprintf('Preprocessing: \n'); 

    if strcmp(dat.modeflag,'classification')
        str{end+1} = sprintf('\t* %s \n', d.PREPROC.groupmode);
    else
        str{end+1} = sprintf('\t* %s \n', d.PREPROC.targetscaling);
    end

    for k=1:numel(d.PREPROC.preprocact)
        str{end+1} = sprintf('\t* Step %g: %s \n', k, d.PREPROC.preprocact{k}); 
    end 

    if FUSION.flag == 3, 
        str = print_modalitydata(str, dat, params, 1); 
    end
    
    if j<nvar, str{end+1} = sprintf('%s \n',repmat('-',1,mxl)); end
    
end

if FUSION.flag ~= 3, str = print_modalitydata(str, dat, params, FUSION.M); end
if nF > 1, str{end+1} = sprintf('%s \n',repmat('=',1,100)); end
str{end+1} = sprintf('Cross-Validation: \n\t* %s\n\n', e.cv);


function str = print_modalitydata(str, dat, params, varind)

e = nk_GetParamDescription2(dat, params.TrainParam.RFE,'FeatFlt');
e = nk_GetParamDescription2(dat, params.TrainParam.RFE,'FeatWrap',e);
e = nk_GetParamDescription2(dat, params.TrainParam,'multiclass',e);
e = nk_GetParamDescription2(dat, params.TrainParam,'GridParam',e);
e = nk_GetParamDescription2(dat, params.TrainParam,'ParamComb',e, varind);
e = nk_GetParamDescription2(dat, params.TrainParam,'SVMprog',e);
e = nk_GetParamDescription2(dat, params.TrainParam,'classifier',e);
e = nk_GetParamDescription2(dat, params.TrainParam,'kernel',e);

% Generate analysis description
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if params.TrainParam.RFE.Filter.flag, 
    str{end+1} = sprintf('Feature selection (Filter): \n\t* %s\n',e.FilterMode);
    str{end+1} = sprintf('\t* %s\n', e.FilterMethod); 
end
if params.TrainParam.RFE.Wrapper.flag,
    str{end+1} = sprintf('Feature selection (Wrapper): \n\t* %s\n',e.WrapperStr);    
    str{end+1} = sprintf('\t* %s\n', e.WrapperMethod); 
end
if params.TrainParam.MULTI.flag, 
    str{end+1} = sprintf('Multi-group classification:'); 
    str{end+1} = sprintf('\n\t* %s\n', e.multiclass);
end
str{end+1} = sprintf('Machine Learning Method: \n\t* %s, %s, %s\n', e.prog, e.classifier, e.kernel);
if e.preML_nCombs>0
    str{end+1} = sprintf('Preprocessing Optimization: %g parameter combinations', e.preML_nCombs); 
    for i=1:numel(e.preML), 
        str{end+1} = sprintf('\n\t(%g) %s', i, e.preML{i}); 
    end; 
    str{end+1} = sprintf('\n');
else
    str{end+1} = sprintf('No optimization of preprocessing parameters\n');
end
if e.ML_nCombs>0
    str{end+1} = sprintf('ML Optimization: %g parameter combinations',e.ML_nCombs); 
    for i=1:numel(e.ML), 
        str{end+1} = sprintf('\n\t(%g) %s',i, e.ML{i}); 
    end; 
    str{end+1} = sprintf('\n');
else
    str{end+1} = sprintf('No optimization of ML parameters\n');
end