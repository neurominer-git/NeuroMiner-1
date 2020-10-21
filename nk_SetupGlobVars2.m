function [status, paramstr] = nk_SetupGlobVars2(dat, act, dispflag, varind)

global ...
    PREPROC ...
    PARMODE ...
    SVM ...
    RFE ...
    CMDSTR ...
    MULTI ...
    SAV ...
    GRD ...
    MODEFL ...
    SCALE ...
    CV ...
    RVM ...
    RAND ...
    VIS ...
    MKLRVM ...
    DATID ...
    SPM5VER ...
    TRAINFUNC ...
    PREDICTFUNC ...
    EVALFUNC ...
    OOCV ...
    LIBSVMTRAIN ...
    LIBSVMPREDICT ...
    FUSION ...
    MULTILABEL ...
    VERBOSE ...
    NM ...
    META ...
    TIME ...
    CVPOS ...
    TEMPL ...
    STACKING

paramstr = [];

switch act
    
    case 'setup_main'
        
        status = 0;
        try 
            CV = dat.cv;
        catch
            paramstr = 'Cross-validation structure';    
        end
        
        try 
            RAND    = dat.TrainParam.RAND;
        catch
            paramstr = sprintf('%s\n%s',paramstr,'Cross-validation setup');
        end
        
        try
            MODEFL  = dat.modeflag;
        catch
            paramstr = sprintf('%s\n%s',paramstr,'Type of predictor: Classification / Regression model');
        end

        try
            FUSION   = dat.TrainParam.FUSION;
        catch
            paramstr = sprintf('%s\n%s',paramstr,'FUSION parameters'); 
        end
        
        try
            STACKING = dat.TrainParam.STACKING;
        catch
            paramstr = sprintf('%s\n%s',paramstr,'STACKING parameters'); 
        end
        
        try
            RAND    = dat.TrainParam.RAND;
        catch
            paramstr = sprintf('%s\n%s',paramstr,'Cross-validation parameters');
        end
        
        try
            SAV     = dat.TrainParam.SAV;
        catch
            paramstr = sprintf('%s\n%s',paramstr,'Save model settings');
        end
 
        if isfield(NM,'OOCV')
            try 
                OOCV = NM.TrainParam.OOCV;
            catch
                paramstr = sprintf('%s\n%s',paramstr,'Independent test validation parameters');
            end
        end
        
        if isempty(NM)
            NM = evalin('base','NM');
        end
        
        if isfield(NM,'time') && ~isempty(NM.time),
            TIME = NM.time;
        end
        
        dat = NM;
        
        if size(dat.label,2)>1
            MULTILABEL.flag = true;
            MULTILABEL.dim  = size(dat.label,2);
        else
            MULTILABEL.flag = false;
            MULTILABEL.dim = 1;
        end
        
        VERBOSE = dat.TrainParam.verbosity;
        
        DATID = dat.id;
        SPM5VER = nk_CheckSPMver;
        
        PARMODE = 0;
        
        if exist('dispflag','var') && ~isempty(dispflag)
            if ~isempty(paramstr), msgbox(paramstr,'Missing parameters detected!','error'); end
        end

    case 'setup_strat'
        
        status = 0;
        if exist('varind','var') && ~isempty(varind)
            
%             if numel(varind)>1
%                 error('CRITICAL: Only one modality can be initialized at a time!!!')
%             elseif isfield(dat,'Y') && varind >= numel(dat.Y);
%                 error('CRITICAL: Modality index exceeds available modalities !!!')
%             end
            
            [PREPROC, ...
                RFE, ...
                GRD, ...
                SVM, ...
                LIBSVMTRAIN, ...
                LIBSVMPREDICT, ...
                RVM, ...
                MKLRVM, ...
                CMDSTR, ...
                MULTI, ...
                VIS, paramstr] = nk_CompatParams2(dat.TrainParam, varind, paramstr);

            if isempty(RFE),        paramstr{end+1} = 'Feature selection parameters'; end
            if isempty(PREPROC),    paramstr{end+1} = 'Preprocessing parameters'; end       
            if isempty(GRD),        paramstr{end+1} = 'Grid optimization settings'; end
            if isempty(MULTI),      paramstr{end+1} = 'Multi-group parameters'; end
            if isempty(VIS),        paramstr{end+1} = 'Visualization parameters'; end
            
            if iscell(PREPROC)
                tPREPROC = PREPROC{1};
            else
                tPREPROC = PREPROC;
            end
            
            if isfield(tPREPROC,'LABELMOD'),
                SCALE.LABELMOD = tPREPROC.LABELMOD; 
            else
                SCALE = [];
            end   
            [TRAINFUNC, PREDICTFUNC] = nk_DefineTrainPredictFunc(true);
            EVALFUNC = nk_DefineEvalFunc;
            
            if isfield(dat.TrainParam,'META')
                META = dat.TrainParam.META;
            else
                META = [];
            end
            
        end
        
        if exist('dispflag','var') && ~isempty(dispflag)
            if ~isempty(paramstr), msgbox(paramstr,'Missing parameters detected!','error'); end
        end
        
    case 'check'
         
        nvar = numel(dat.Y); paramstr = [];
        status = 0;
            
        checkfields = {'PREPROC', 'SVM',  'GRD',  'RFE', 'VIS'};
        descriptors = {'Preprocessing', ...
                        'ML algorithm', ...
                        'ML optimization', ...
                        'Feature selection', ...
                        'Visualization'};
        if numel(unique(NM.label))>2 && strcmp(NM.modeflag,'classification')
            checkfields = [checkfields 'MULTI'];
            descriptors = [descriptors 'Multi-group'];
        end
        
        if isfield(dat,'TrainParam') && ~isempty(dat.TrainParam)
             if  nvar > 1, 
                 if ~isfield(dat.TrainParam,'FUSION') || isempty(dat.TrainParam.FUSION)
                    paramstr = 'Fusion parameters'; status = 1;
                 else
                     switch dat.TrainParam.FUSION.flag 
                         case {0,1,2}
                             checkfields = [checkfields, 'SAV', 'RAND']; 
                             descriptors = [descriptors, 'Saving parameters', 'Cross-validation parameters'];
                             params = dat.TrainParam; 
                            [status, paramstr] = check_params(paramstr, params, checkfields, descriptors);
                            if status
                                paramstr = sprintf('Missing Parameters: %s', paramstr(3:end));
                            end
                         case 3
                             paramstr = cell(numel(dat.TrainParam.FUSION.M),1); ll = 1; 
                             for i = 1:nvar
                                 if ~sum(any(dat.TrainParam.FUSION.M == i)), continue; end
                                 params = dat.TrainParam.STRAT{i};
                                 [lstatus, paramstr{ll}] = check_params(paramstr{ll}, params, checkfields, descriptors);
                                 if lstatus
                                    status = 1; paramstr{ll} = sprintf('Modality #%g: Missing parameters =>%s',i, paramstr{ll}(2:end));
                                 else
                                    paramstr{ll} = sprintf('Modality #%g: OK',i); 
                                 end
                                 ll=ll+1;
                             end
                            
                     end
                     
                 end
             else
                 checkfields = [checkfields, 'SAV', 'RAND'];
                 descriptors = [descriptors, 'Saving parameters', 'Cross-validation definitions'];
                 params = dat.TrainParam; 
                 [status, paramstr] = check_params(paramstr, params, checkfields, descriptors);
             end     
        else
            paramstr = 'Training parameters'; status = 1;
        end
        if ~status, 
            paramstr = [];
        elseif iscell(paramstr), 
            paramstr = char(paramstr); 
        end
        if ~isfield(dat,'modeflag')
           paramstr = char(paramstr,'Prediction framework undefined'); status = 1;
        end
        if ~isfield(dat,'cv') 
            paramstr = char(paramstr,'Cross-validation structure undefined');  status = 1;  
        end 
        if exist('dispflag','var') && ~isempty(dispflag)
            if ~isempty(paramstr) && status
                msgbox(paramstr,'INCOMPLETE SETUP!','error');
            end
        end
        
    case 'clear'
        clear   global ...
                PREPROC ...
                PARMODE ...
                SVM ...
                RFE ...
                CMDSTR ...
                MULTI ...
                TEST ...
                SAV ...
                GRD ...
                MODEFL ...
                CV ...
                RVM ...
                RAND ...
                VIS ...
                MKLRVM ...
                DATID ...
                SPM5VER ...
                TRAINFUNC ...
                PREDICTFUNC ...
                EVALFUNC ...
                OOCV ...
                LIBSVMTRAIN ...
                LIBSVMPREDICT ...
                VERBOSE ...
                TEMPL ...
                META ...
                STACKING ...
                CVPOS ...
                TIME
    otherwise
        error(['Option ' act ' not available!'])
end
end

function [status, paramstr] = check_params(paramstr, params, checkfields, descriptors)

istr = []; status = 0;
for i=1:numel(checkfields)
    if ~isfield(params,checkfields{i}) || isempty(params.(checkfields{i}))
        istr = sprintf('%s, %s', istr, descriptors{i}); status = 1; 
    end
end
if ~isempty(istr)
    if ~isempty(paramstr)
        paramstr = [paramstr istr];
    else
        paramstr = istr;
    end
end
end