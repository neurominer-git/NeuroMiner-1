function nk_PrintWs(dat, oTrainParam)

clc
if iscell(dat.Y), nvar = length(dat.Y); else, nvar = 1; end
showmodalmax = 3;
showmodalvec = 1:nvar;
nF=1;
if isempty(oTrainParam.FUSION), oTrainParam.FUSION.flag = 0; end
if nvar > 1
    if oTrainParam.STACKING.flag == 1
        cprintf('*black', 'WORKSPACE VIEWER: ');
        cprintf('*black', 'STACKING OPERATES ON ANALYSES %s ', strjoin(string(oTrainParam.STACKING.sel_anal),', ')); 
    else
        nF = numel(oTrainParam.FUSION.M);
        if nF>1,
            strmod = sprintf('%g MODALITIES:',nF);
        else
            strmod = 'MODALITY';
        end
        cprintf('*black', 'WORKSPACE VIEWER: %s ', strmod); 
        for j=1:nF, 
            if j > showmodalmax, fprintf(' ...'); break; end;
            cprintf('*black','#%g ', oTrainParam.FUSION.M(j)); 
        end
    end
    fprintf('\n') 
end

 e = nk_GetParamDescription2(dat, dat,'cv');

% Loop through variates
for jj=1:nF

    j = showmodalvec(jj);

    if oTrainParam.FUSION.flag == 1 && j>1
        continue; 
    elseif oTrainParam.FUSION.flag == 3
        TrainParam = oTrainParam.STRAT{j};
    else
        TrainParam = oTrainParam;
    end

    switch oTrainParam.FUSION.flag
        case {0,2}
            % Get info about modality j
            d = nk_GetParamDescription2(dat, dat, 'VarDesc', [], oTrainParam.FUSION.M(j));
            d = nk_GetParamDescription2(dat, TrainParam.PREPROC{oTrainParam.FUSION.M(j)},'PreProc', d, j);
        case 1
            % Get info about all modalities
            d = nk_GetParamDescription2(dat, dat,'VarDesc', [], oTrainParam.FUSION.M);
            d = nk_GetParamDescription2(dat, TrainParam.PREPROC{j},'PreProc', d, 1);
        case 3
            % Get info about modality j
            d = nk_GetParamDescription2(dat, dat,'VarDesc', [], oTrainParam.FUSION.M);
            d = nk_GetParamDescription2(dat, TrainParam.PREPROC,'PreProc', d);                    
    end
    mxl = 0;
    if oTrainParam.STACKING.flag == 1
        mdlstr = sprintf('ANALYSIS SETUP FOR STACKING\n');

    else
        if oTrainParam.FUSION.flag == 1
            for i=1:numel(d.datadescriptor)
                mdlstr = sprintf('MODALITY %g : %s \n', oTrainParam.FUSION.M(i), d.datadescriptor{i});
                cprintf('*blue',mdlstr);
                mxli = size(mdlstr,2);
                if mxli > mxl, mxl = mxli; end
            end
        else
            if numel(d.datadescriptor) > 1,
                datdesc = d.datadescriptor{j};
            else
                datdesc = d.datadescriptor{1};
            end
            mdlstr = sprintf('MODALITY %g : %s \n', oTrainParam.FUSION.M(j), datdesc);
        end
    end
    cprintf('*blue',mdlstr);
    mxl = size(mdlstr,2);
    cprintf('*blue','%s \n',repmat('*',1,mxl));
    if oTrainParam.STACKING.flag == 1
        cprintf('*black','Input layer analyses: ')
        for i=1:numel(oTrainParam.STACKING.featnames)
            fprintf('\n\t* %s',oTrainParam.STACKING.featnames{i})
        end
        fprintf('\n');
        predstr = {'CV1 training data','CV1 test data'};
        cprintf('*black','Prediction extraction for stacking: ')
        fprintf('\n\t* from %s\n', predstr{oTrainParam.STACKING.mode});
    end
    cprintf('*black','Preprocessing: \n'); 

    if strcmp(dat.modeflag,'classification')
        fprintf('\t* %s \n', d.PREPROC.groupmode);
    else
        fprintf('\t* %s \n', d.PREPROC.targetscaling);
    end

    for k=1:numel(d.PREPROC.preprocact)
        fprintf('\t* Step %g: %s \n', k, d.PREPROC.preprocact{k}); 
    end 

    if oTrainParam.FUSION.flag == 3, print_modalitydata(dat, TrainParam, 1); end

    if jj > showmodalmax, 
        cprintf('*blue','%s \n',repmat('-',1,mxl)); 
        cprintf('*red','>>> %g further modalities included in this analysis... \n', nF-j+1);
        break; 
    else
        if oTrainParam.STACKING.flag ~= 1
            if j<nvar, cprintf('*blue','%s \n',repmat('-',1,mxl)); end
        end
    end
end
if oTrainParam.FUSION.flag ~= 3, print_modalitydata(dat, oTrainParam, oTrainParam.FUSION.M); end
if nF > 1, cprintf('*blue','%s \n',repmat('=',1,100)); end
cprintf('*black','Cross-Validation: '); 
fprintf('\n\t* %s\n\n', e.cv);
fprintf('Press any key to return to NM ...')
pause

function print_modalitydata(dat, TrainParam, varind)

e = nk_GetParamDescription2(dat, TrainParam.RFE,'FeatFlt');
e = nk_GetParamDescription2(dat, TrainParam.RFE,'FeatWrap',e);
e = nk_GetParamDescription2(dat, TrainParam,'multiclass',e);
e = nk_GetParamDescription2(dat, TrainParam,'GridParam',e);
e = nk_GetParamDescription2(dat, TrainParam,'ParamComb',e, varind);
e = nk_GetParamDescription2(dat, TrainParam,'SVMprog',e);
e = nk_GetParamDescription2(dat, TrainParam,'classifier',e);
e = nk_GetParamDescription2(dat, TrainParam,'kernel',e);

% Generate analysis description
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if TrainParam.RFE.Filter.flag, 
    cprintf('*black','Feature selection (Filter): '); fprintf('\n\t* %s\n',e.FilterMode);
    fprintf('\t* %s\n', e.FilterMethod); 
end
if TrainParam.RFE.Wrapper.flag,
    cprintf('*black','Feature selection (Wrapper): '); fprintf('\n\t* %s\n',e.WrapperStr);    
    fprintf('\t* %s\n', e.WrapperMethod); 
end
if TrainParam.MULTI.flag, 
    cprintf('*black','Multi-group classification: '); fprintf('\n\t* %s\n', e.multiclass);
end
cprintf('*black','Machine Learning Method: '); fprintf('\n\t* %s, %s, %s\n', e.prog, e.classifier, e.kernel);
if e.preML_nCombs>0
    cprintf('*black','Preprocessing Optimization: %g parameter combinations', e.preML_nCombs); 
    for i=1:numel(e.preML), fprintf('\n\t(%g) %s', i, e.preML{i}); end; fprintf('\n');
else
    cprintf('*black','No optimization of preprocessing parameters \n')
end
if e.ML_nCombs>0
    cprintf('*black','ML Optimization: %g parameter combinations',e.ML_nCombs); 
    for i=1:numel(e.ML), fprintf('\n\t(%g) %s',i, e.ML{i}); end; fprintf('\n');
else
    cprintf('*black','No optimization of ML parameters \n')
end