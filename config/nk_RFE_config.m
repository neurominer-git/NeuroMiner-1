function [act, RFE] = nk_RFE_config(act, TrainParam, SVM, MODEFL, MULTI, GRD, defaultsfl, parentstr)

if ~exist('defaultsfl','var') || isempty(defaultsfl),  defaultsfl = 0; end

if ~isempty(act) || ~defaultsfl

    if ~isfield(TrainParam,'RFE') || isempty(TrainParam.RFE); 
        [~, TrainParam.RFE ] = nk_RFE_config([], TrainParam, SVM, MODEFL, MULTI, true);
    end
    if ~isfield(TrainParam.RFE,'Filter')
        [~, TrainParam.RFE ] = nk_Filter_config([], TrainParam.RFE, SVM, MODEFL, MULTI, 1); 
    end
    if ~isfield(TrainParam.RFE,'Wrapper')
        [~, TrainParam.RFE ] = nk_Wrapper_config([], TrainParam.RFE, SVM, MODEFL, 1); 
    end
    
    RFE = TrainParam.RFE;
    RFE.dispres = 1;
    nk_PrintLogo
    
    immretstr = ''; menuvec = 1:3;
   
    if prod(GRD.n_params) <=1 && ...
       ( isfield(TrainParam,'RFE') && isfield(TrainParam.RFE,'Wrapper')   && ~TrainParam.RFE.Wrapper.flag ) && ...
       ( isfield(TrainParam,'RFE') && isfield(TrainParam.RFE,'Filter')    && ~TrainParam.RFE.Filter.flag )
        if isfield(RFE,'RetrainImmediate') && RFE.RetrainImmediate, immstr = 'yes'; else immstr = 'no'; end
        immretstr = sprintf('Skip CV1 cycle and train on full CV1 training and test data [ %s ]|',immstr);
        menuvec = 1:4;
    else
        % Deactivate immediate retrain!
        if isfield(RFE,'RetrainImmediate') && RFE.RetrainImmediate, cprintf('red','Full CV1 partition training DISABLED due to optimization requirements'); end
        RFE.RetrainImmediate = 0;
    end
    
    STATUS = nk_CheckFieldStatus(TrainParam,{'RFE'},{'Filter','Wrapper','CV2Class'});
    mestr = 'Ensemble generation setup'; navistr = sprintf('%s\n\t>>> %s',parentstr, mestr); cprintf('*blue','\nYou are here: %s >>> ',parentstr);
    
    act = nk_input(mestr,0, 'mq', ...
                    [['Filter-based ensemble generation [ ' STATUS.Filter ' ]|'] ...
                    ['Wrapper-based ensemble generation [ ' STATUS.Wrapper ' ]|'] ...
                    ['More ensemble learning options [ ' STATUS.CV2Class ' ]|'] ...
                    immretstr ],menuvec);
    switch act
        case 1
            if ~isfield(RFE,'Filter') || isempty(RFE.Filter), [~, RFE] = nk_Filter_config([], RFE, SVM, MODEFL, MULTI, 1); end
            acti = 1; while acti>0, [acti, RFE] = nk_Filter_config(acti, RFE, SVM, MODEFL, MULTI, [], navistr); end
        case 2
            if ~isfield(RFE,'Wrapper') || isempty(RFE.Wrapper), [~, RFE] = nk_Wrapper_config([], RFE, SVM, MODEFL, GRD, 1); end
            acti = 1; while acti>0, [acti, RFE] = nk_Wrapper_config(acti, RFE, SVM, MODEFL, GRD, MULTI, [], navistr); end
        case 3
            RFE = config_ens(RFE, MODEFL, navistr);    
        case 4
            RFE.RetrainImmediate = nk_input('Train models using entire CV1 data partition?',0,'yes|no',[1,0],1);
    end
else
    RFE = [];
    [~,RFE] = nk_Filter_config([], RFE, SVM, MODEFL, MULTI,1);
    [~,RFE] = nk_Wrapper_config([], RFE, SVM, MODEFL, GRD, MULTI, 1);
    RFE.CV2Class.type = 2;
    RFE.CV2Class.EnsembleStrategy.type = 0;
    RFE.CV2Class.EnsembleStrategy.DataType = 2;
    RFE.CV2Class.EnsembleStrategy.Metric = 2;
    RFE.CV2Class.EnsembleStrategy.AggregationLevel = 2;
    RFE.ClassRetrain = 1;
    RFE.RetrainImmediate = 0;
    RFE.dispres = 1;
    act =0;
end

TrainParam.RFE = RFE;

end

% *************************** CV2 Classification **************************
function RFE = config_ens(RFE, MODEFL, parentstr)
% This submodule configures how the class memberships or the continuous labels 
% of the CV2 test data should be predicted.
        
% Use always the ensemble prediction
RFE.CV2Class.type = 2;

flgWrp = 0; flgFlt = 0;
if isfield(RFE,'Filter') && isfield(RFE.Filter,'EnsembleStrategy') && ...
        isfield(RFE.Filter.EnsembleStrategy,'type') && ...
        RFE.Filter.EnsembleStrategy.type ~= 9
    flgFlt = 1;
end

if isfield(RFE,'Wrapper') && isfield(RFE.Wrapper,'EnsembleStrategy') &&...
        isfield(RFE.Wrapper.EnsembleStrategy,'type') && ...
        RFE.Wrapper.EnsembleStrategy.type ~= 9
    flgWrp = 1;
end

CV2Mode = 1;
% if flgWrp || flgFlt
%     CV2Mode = nk_input('Define CV2 ensemble construction method',0, 'm',...
%                          ['Aggregate all CV1 ensembles of given CV2 partition|' ...
%                           'Take over optimization parameters of CV1 ensemble construction (not functional)|' ...
%                           'Define independent CV2 ensemble construction strategy (not functional)'],1:3);
% else
%     CV2Mode = nk_input('Define CV2 ensemble construction method',0, 'm',...
%                          ['Aggregate all CV1 ensembles of given CV2 partition|' ...
%                           'Define independent CV2 ensemble construction strategy (not functional)'],[1,3]);
%     
% end

switch CV2Mode
   
    case 1
        
        RFE.CV2Class.EnsembleStrategy.type = 0;
        RFE.CV2Class.EnsembleStrategy.DataType = 0;
        
        if flgWrp || flgFlt
            if flgWrp, 
                RFE.CV2Class.EnsembleStrategy.Metric = RFE.Wrapper.EnsembleStrategy.Metric;
            else
                RFE.CV2Class.EnsembleStrategy.Metric = RFE.Filter.EnsembleStrategy.Metric;
            end
        else
            switch MODEFL
                case 'classification'
                    RFE.CV2Class.EnsembleStrategy.Metric = spm_input('Metric',0,'m', ...
                      'Target matrix (majority voting)||Decision (SVM) / Probability (RVM) value matrix',1:2); 
                case 'regression'
                    RFE.CV2Class.EnsembleStrategy.Metric = 2;
            end
        end    
        
    case 2
        if flgWrp
            RFE.CV2Class.EnsembleStrategy = RFE.Wrapper.EnsembleStrategy;
        elseif flgFlt
            RFE.CV2Class.EnsembleStrategy = RFE.Filter.EnsembleStrategy;
        end
    case 3
        RFE.CV2Class = nk_EnsembleStrategy2_config(NM, RFE.CV2Class, [], parentstr);
end
    
% Should the base learners be retrained using the entire CV1 data ? 
RFE.ClassRetrain = ...
    nk_input('Retrain classifiers with all data (CV1 training + CV1 test data) ?',0,'yes|no',[1,0],1);

% Should the CV1 ensembles' outputs be aggregated or their base learners
% outputs be pulled together across the CV2 partitions in order to create 
% an CV2-OOT ensemble
RFE.CV2Class.EnsembleStrategy.AggregationLevel = ...
    nk_input('Aggregating predictions across CV2 permutations',0,'m', ...
    ['Mean ensemble prediction across CV2 permutations (grand mean)|' ...
    'Grand ensemble prediction constructed across CV2 permutations'],[0,1],2);

end
