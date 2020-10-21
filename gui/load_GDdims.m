function handles = load_GDdims(handles, Params, Label, GDdims)

handles.MLparams = []; handles.grid = [];

% Add binary-classification data to handles
if isfield( GDdims,'Model'), 
    handles.MLparams = GDdims.Model;
elseif isfield( GDdims,'MLparams' )
    handles.MLparams = GDdims.MLparams;
else
    handles.MLparams.NumParamCombs = 1;
end
if isfield(GDdims,'grid')
    handles.grid     = GDdims.grid;    
end

if strcmp(handles.selYaxis.String{handles.selYaxis.Value},'Multi-group probabilities derived from similarity averaging')
    fld = 'MultiClassProb';
else
    fld = 'MultiClass';
end

if isfield(handles,'SubIndex') && ~handles.oocvview, I = handles.SubIndex; else, I = true(size(handles.NM.label,1),1); end

if strcmp(handles.METAstr,'none')
    handles.ModelParams                         = GDdims.Model.ParamCombs;
    handles.ModelParamsDesc                     = GDdims.Model.ParamDesc;
end

% Check type of analysis
if isfield(GDdims,'BinClass') || isfield(GDdims,'MultiClass') 
    
    handles.modeflag                            = 'classification';
    handles.nclass                              = numel(GDdims.BinClass);
    
    for j=1:handles.nclass
        handles.BinClass{j}.description         = Params.class{1,1}{j}.groupdesc; %NM.cv.class
        handles.BinClass{j}.groupind            = Params.class{1,1}{j}.groups;
        handles.BinClass{j}.groupnames          = regexp(handles.BinClass{j}.description,' vs ','split');
    end
    [handles.labels, handles.sortind ]          = sort(Label(:,handles.curlabel),'ascend');

    for j=1:handles.nclass

        handles.BinClass{j}.ind                 = GDdims.BinClass{j}.index_predictions & I;
        handles.BinClass{j}.cases               = handles.subjects(handles.BinClass{j}.ind);
        switch handles.tglSort.Value
            case 1
                switch handles.selYaxis.String{handles.selYaxis.Value}
                    case {'Mean classifier scores','Mean classifier scores (95%-CIs)','Mean classifier scores (SD)'}
                        P = GDdims.BinClass{j}.mean_predictions( handles.BinClass{j}.ind );
                    case 'Ensemble-based probability scores'
                        P = GDdims.BinClass{j}.prob_predictions( handles.BinClass{j}.ind ,1 );
                    case {'Mean ensemble-based probability scores (95%-CIs)','Mean ensemble-based probability scores (SD)'}
                        P = GDdims.CV2grid.mean_predictions( handles.BinClass{j}.ind , j, handles.curlabel);
                end
                [ ~, ...
                        handles.BinClass{j}.sortind] = sort(P,'ascend');
                handles.BinClass{j}.labels = Label(handles.BinClass{j}.ind, handles.curlabel);
                handles.BinClass{j}.labels = handles.BinClass{j}.labels(handles.BinClass{j}.sortind);
               
            case 0
                [ handles.BinClass{j}.labels, ...
                    handles.BinClass{j}.sortind]        = sort(Label(handles.BinClass{j}.ind, handles.curlabel),'ascend');
        end
        handles.BinClass{j}.cases               = handles.BinClass{j}.cases(handles.BinClass{j}.sortind);
        handles.BinClass{j}.labelh              = zeros(size(handles.BinClass{j}.labels,1),1);
        if numel(handles.BinClass{j}.groupind) == 2
            handles.BinClass{j}.ind1 = handles.BinClass{j}.labels == handles.BinClass{j}.groupind(1); 
            handles.BinClass{j}.ind2 = handles.BinClass{j}.labels == handles.BinClass{j}.groupind(2);
            handles.BinClass{j}.labelh(handles.BinClass{j}.ind1,1) = 1; 
            handles.BinClass{j}.labelh(handles.BinClass{j}.ind2,1) = -1;
            handles.BinClass{j}.one_vs_all = false;
        else
            handles.BinClass{j}.ind1 = handles.BinClass{j}.labels == handles.BinClass{j}.groupind(1); 
            handles.BinClass{j}.ind2 = handles.BinClass{j}.labels ~= handles.BinClass{j}.groupind(1);
            handles.BinClass{j}.groupind(2) = handles.nclass;
            handles.BinClass{j}.labelh(handles.BinClass{j}.ind1,1) = 1; 
            handles.BinClass{j}.labelh(~handles.BinClass{j}.ind1,1) = -1;
            handles.BinClass{j}.one_vs_all = true;
        end
     
        if isfield(GDdims,'CV2grid') && isfield(GDdims,'mean_predictions')
            handles.BinClass{j}.CV2grid.mean_predictions    = GDdims.CV2grid.mean_predictions( handles.BinClass{j}.ind, j, handles.curlabel );
            handles.BinClass{j}.CV2grid.std_predictions     = GDdims.CV2grid.std_predictions( handles.BinClass{j}.ind, j, handles.curlabel  );
            handles.BinClass{j}.CV2grid.CI1_predictions     = GDdims.CV2grid.CI_predictions( handles.BinClass{j}.ind, 1, j, handles.curlabel );
            handles.BinClass{j}.CV2grid.CI2_predictions     = GDdims.CV2grid.CI_predictions( handles.BinClass{j}.ind, 2, j, handles.curlabel );
            handles.BinClass{j}.CV2grid.mean_predictions    = handles.BinClass{j}.CV2grid.mean_predictions ( handles.BinClass{j}.sortind );
            handles.BinClass{j}.CV2grid.std_predictions     = handles.BinClass{j}.CV2grid.std_predictions ( handles.BinClass{j}.sortind );
            handles.BinClass{j}.CV2grid.CI1_predictions     = handles.BinClass{j}.CV2grid.CI1_predictions ( handles.BinClass{j}.sortind );
            handles.BinClass{j}.CV2grid.CI2_predictions     = handles.BinClass{j}.CV2grid.CI2_predictions ( handles.BinClass{j}.sortind );
        end
        handles.BinClass{j}.mean_predictions    = GDdims.BinClass{j}.mean_predictions( handles.BinClass{j}.ind );
        handles.BinClass{j}.std_predictions     = GDdims.BinClass{j}.std_predictions( handles.BinClass{j}.ind );
        handles.BinClass{j}.CI1_predictions     = GDdims.BinClass{j}.CI1_predictions( handles.BinClass{j}.ind );
        handles.BinClass{j}.CI2_predictions     = GDdims.BinClass{j}.CI2_predictions( handles.BinClass{j}.ind  );
        handles.BinClass{j}.prob_predictions    = GDdims.BinClass{j}.prob_predictions( handles.BinClass{j}.ind ,:);
        if isfield(GDdims.BinClass{j},'mean_cutoff_probabilities')
            handles.BinClass{j}.mean_cutoff_probabilities = GDdims.BinClass{j}.mean_cutoff_probabilities( handles.BinClass{j}.ind );
            handles.BinClass{j}.std_cutoff_probabilities =  GDdims.BinClass{j}.std_cutoff_probabilities( handles.BinClass{j}.ind );
            handles.BinClass{j}.mean_cutoff_percentiles = GDdims.BinClass{j}.mean_cutoff_percentiles( handles.BinClass{j}.ind );
            handles.BinClass{j}.std_cutoff_percentiles =  GDdims.BinClass{j}.std_cutoff_percentiles( handles.BinClass{j}.ind );
            handles.BinClass{j}.mean_cutoff_probabilities  = handles.BinClass{j}.mean_cutoff_probabilities( handles.BinClass{j}.sortind, handles.curlabel );
            handles.BinClass{j}.std_cutoff_probabilities = handles.BinClass{j}.std_cutoff_probabilities ( handles.BinClass{j}.sortind, handles.curlabel );
            handles.BinClass{j}.mean_cutoff_percentiles = handles.BinClass{j}.mean_cutoff_percentiles ( handles.BinClass{j}.sortind, handles.curlabel );
            handles.BinClass{j}.std_cutoff_percentiles = handles.BinClass{j}.std_cutoff_percentiles ( handles.BinClass{j}.sortind, handles.curlabel );
            handles.BinClass{j}.mean_globalcutoff_probabilities = GDdims.BinClass{j}.mean_globalcutoff_probabilities ;
            handles.BinClass{j}.std_globalcutoff_probabilities = GDdims.BinClass{j}.std_globalcutoff_probabilities ;
            handles.BinClass{j}.mean_globalcutoff_percentiles = GDdims.BinClass{j}.mean_globalcutoff_percentiles ;
            handles.BinClass{j}.std_globalcutoff_percentiles= GDdims.BinClass{j}.std_globalcutoff_percentiles;
            if isfield(GDdims.BinClass{j},'globalcutoff_probabilities')
                handles.BinClass{j}.globalcutoff_probabilities= GDdims.BinClass{j}.globalcutoff_probabilities( handles.BinClass{j}.ind );
            end
            handles.BinClass{j}.CoxMode = 1;
        else
            handles.BinClass{j}.CoxMode = 0;
        end
        handles.BinClass{j}.mean_predictions    = handles.BinClass{j}.mean_predictions( handles.BinClass{j}.sortind, handles.curlabel );
        handles.BinClass{j}.std_predictions     = handles.BinClass{j}.std_predictions( handles.BinClass{j}.sortind, handles.curlabel );
        handles.BinClass{j}.CI1_predictions     = handles.BinClass{j}.CI1_predictions( handles.BinClass{j}.sortind, handles.curlabel );
        handles.BinClass{j}.CI2_predictions     = handles.BinClass{j}.CI2_predictions( handles.BinClass{j}.sortind, handles.curlabel );
        handles.BinClass{j}.prob_predictions    = handles.BinClass{j}.prob_predictions( handles.BinClass{j}.sortind ,:);
        handles.BinClass{j}.prob_finalpred      = GDdims.BinClass{j}.prob_finalpred( handles.BinClass{j}.sortind, handles.curlabel );
      
        if handles.BinClass{j}.CoxMode
            if isfield(GDdims.BinClass{j},'globalcutoff_probabilities')
                handles.BinClass{j}.globalcutoff_probabilities = handles.BinClass{j}.globalcutoff_probabilities( handles.BinClass{j}.sortind ,:);
                handles.BinClass{j}.prob_contingency = ALLPARAM(handles.BinClass{j}.labelh, handles.BinClass{j}.prob_predictions(:,1)-handles.BinClass{j}.globalcutoff_probabilities);
            else
                handles.BinClass{j}.prob_contingency = GDdims.BinClass{j}.contigency;
            end
            handles.BinClass{j}.contingency      = handles.BinClass{j}.prob_contingency;
        else
            handles.BinClass{j}.contingency         = ALLPARAM(handles.BinClass{j}.labelh, handles.BinClass{j}.mean_predictions);
            handles.BinClass{j}.prob_contingency    = ALLPARAM(handles.BinClass{j}.labelh, handles.BinClass{j}.prob_predictions(:,1)-0.5);
        end
        switch handles.METAstr
            case 'none'
                handles.BinClass{j}.best_TR             = GDdims.bestTR{j}(:,:, handles.curlabel);
                handles.BinClass{j}.best_TS             = GDdims.bestTS{j}(:,:, handles.curlabel);
                handles.BinClass{j}.best_P              = GDdims.bestP{j}(:,:, handles.curlabel);
                handles.BinClass{j}.best_Ppos           = GDdims.bestPpos{j}(:, handles.curlabel);
                handles.BinClass{j}.best_CVperf         = GDdims.best_CVperf(j);
                handles.BinClass{j}.best_TSperf         = GDdims.best_TSperf(j);
                handles.BinClass{j}.best_Complexity     = GDdims.bestComplexity{j};
                handles.BinClass{j}.best_Error          = GDdims.bestError{j};
        end
        
        % For ROC analysis
        [handles.BinClass{j}.X , ...
         handles.BinClass{j}.Y ] = perfcurve2(handles.BinClass{j}.labelh, handles.BinClass{j}.mean_predictions, 1);
        
        % Prepare table for export
        handles.BinClass{j}.tbl.colnames        = {'Cases', 'EXP_LABEL', 'PRED_LABEL', 'Errors', 'Mean_Score', 'Std_Score', 'Ens_ProbPred'};
        handles.BinClass{j}.tbl.rownames        = handles.BinClass{j}.cases;
        handles.BinClass{j}.tbl.array           = [handles.BinClass{j}.labelh ...
                                                   handles.BinClass{j}.prob_finalpred ...
                                                   handles.BinClass{j}.labelh ~= handles.BinClass{j}.prob_finalpred ...
                                                   handles.BinClass{j}.mean_predictions ...
                                                   handles.BinClass{j}.std_predictions ...
                                                   handles.BinClass{j}.prob_predictions(:,1)];
        handles.BinClass{j}.tbl_cont.rownames   = fieldnames(handles.BinClass{j}.contingency);
        handles.BinClass{j}.tbl_cont.colnames   = {'Metric', handles.BinClass{j}.description};
        handles.BinClass{j}.tbl_cont.array      = struct2cell( handles.BinClass{j}.contingency);
        remind = find(strcmp('FPRvec',handles.BinClass{j}.tbl_cont.rownames) | strcmp('TPRvec',handles.BinClass{j}.tbl_cont.rownames) | strcmp('X',handles.BinClass{j}.tbl_cont.rownames));
        handles.BinClass{j}.tbl_cont.array(remind) = [];
        handles.BinClass{j}.tbl_cont.array = cell2mat(handles.BinClass{j}.tbl_cont.array);
        handles.BinClass{j}.tbl_cont.rownames(remind) = [];
        
    end    
    % Add multi-class data to handles
    if isfield(GDdims,'MultiClass')
        handles.ngroups                             = numel(handles.NM.groupnames);
        handles.MultiClass                          = GDdims.(fld);
        indn                                        = ~isnan(GDdims.(fld).multi_probabilitiesCV2(:,1,handles.curlabel)) & I;
        [~, handles.MultiClass.sortind]             = sort(Label(:,handles.curlabel),'ascend');
        indn                                        = handles.MultiClass.sortind(indn(handles.MultiClass.sortind));
        handles.MultiClass.cases                    = handles.subjects(indn);
        handles.MultiClass.labels                   = Label(indn, handles.curlabel);    
        handles.MultiClass.probabilities            = GDdims.(fld).multi_probabilitiesCV2(indn , :, handles.curlabel);
        handles.MultiClass.onevsall_labels          = zeros(numel(indn),handles.ngroups);
        handles.MultiClass.onevsall_scores          = zeros(numel(indn),handles.ngroups);
        for j = 1:handles.ngroups
            ind = true(1,handles.ngroups); ind(j) = false;
            probrest = nanmean(handles.MultiClass.probabilities(:,ind),2);
            %probone  = handles.MultiClass.probabilities(:,j);
            handles.MultiClass.onevsall_labels(:,j) = handles.MultiClass.labels == j;     
            handles.MultiClass.onevsall_scores(:,j) = 1-probrest;
        end
        handles.MultiClass.onevsall_labels(handles.MultiClass.onevsall_labels == 0) = -1;
        for j = 1:handles.ngroups
            [handles.MultiClass.X{j}, ...
                handles.MultiClass.Y{j}, ...
                handles.MultiClass.T{j}, ...
                handles.MultiClass.class{j}.AUC] = ...
                            perfcurve2(handles.MultiClass.onevsall_labels(:,j), handles.MultiClass.onevsall_scores(:,j), 1);
        end
        handles.MultiClass.errors                   = handles.MultiClass.errors(indn);
        handles.MultiClass.multi_predictionsCV2     = GDdims.(fld).multi_predictionsCV2(indn, handles.curlabel);
        handles.MultiClass.std_multi_predictionsCV2 = GDdims.(fld).multi_predictionsCV2_std(indn, handles.curlabel);
        handles.MultiClass.CI1_multi_predictionsCV2 = GDdims.(fld).multi_predictionsCV2_ci1(indn, handles.curlabel);
        handles.MultiClass.CI2_multi_predictionsCV2 = GDdims.(fld).multi_predictionsCV2_ci2(indn, handles.curlabel);
         switch handles.METAstr
            case 'none'
                handles.MultiClass.best_TR                  = GDdims.multi_bestTR;
                handles.MultiClass.best_TS                  = GDdims.multi_bestTS;
                handles.MultiClass.best_CVperf              = GDdims.best_MultiCVperf;
                handles.MultiClass.std_best_CVperf          = GDdims.best_sdMultiCVperf;
                handles.MultiClass.best_TSperf              = GDdims.best_MultiTSperf;
                handles.MultiClass.std_best_TSperf          = GDdims.best_sdMultiTSperf;
                handles.MultiClass.best_Error               = handles.MultiClass.best_TR - handles.MultiClass.best_TS;
         end
        Compl                                       = zeros(size(handles.BinClass{1}.best_Complexity));
        for curclass = 1:numel(handles.BinClass)
            Compl = Compl + handles.BinClass{curclass}.best_Complexity;
        end
        Compl = Compl./numel(handles.BinClass);
        handles.MultiClass.best_Complexity          = Compl;
        % Prepare table for export
        handles.MultiClass.tbl.colnames        = {'Cases', 'EXP_LABEL', 'PRED_LABEL', 'Errors'};
        handles.MultiClass.tbl.rownames        = handles.MultiClass.cases;
        handles.MultiClass.tbl.array           = [handles.MultiClass.labels ...
                                                   handles.MultiClass.multi_predictionsCV2 ....
                                                   handles.MultiClass.errors' ];                                  
        for j = 1:handles.ngroups
            handles.MultiClass.tbl.colnames    = [ handles.MultiClass.tbl.colnames sprintf('Probabilities_%g',j) ];                          
            handles.MultiClass.tbl.array       = [ handles.MultiClass.tbl.array handles.MultiClass.probabilities(:,j) ];                                        
        end
        
        for j = 1:handles.ngroups
            handles.MultiClass.tbl.colnames    = [ handles.MultiClass.tbl.colnames ...
                                                    sprintf('EXP_%s_vs_REST',handles.NM.groupnames{j}) ...
                                                    sprintf('PRED_%s_vs_REST',handles.NM.groupnames{j}) ...
                                                    sprintf('Score_%s_vs_REST',handles.NM.groupnames{j})];  
            handles.MultiClass.tbl.array       = [ handles.MultiClass.tbl.array  ...
                                                    handles.MultiClass.onevsall_labels(:,j) ...
                                                    sign(handles.MultiClass.onevsall_scores(:,j)) ...
                                                    handles.MultiClass.onevsall_scores(:,j) ];  
        end
        handles.MultiClass.tbl_cont.rownames   = fieldnames(handles.MultiClass.class{1});
        handles.MultiClass.tbl_cont.colnames   = {'Metric'};
        handles.MultiClass.tbl_cont.array      = [];
        remind = find(  strcmp('FPRvec',handles.MultiClass.tbl_cont.rownames) | ...
                        strcmp('TPRvec',handles.MultiClass.tbl_cont.rownames) | ...
                        strcmp('X',handles.MultiClass.tbl_cont.rownames));
        for j = 1:handles.ngroups, 
            handles.MultiClass.tbl_cont.colnames{j+1} = sprintf('%s vs REST',handles.NM.groupnames{j});
            arr = struct2cell( handles.MultiClass.class{j});
            arr(remind)=[];
            handles.MultiClass.tbl_cont.array = [handles.MultiClass.tbl_cont.array arr];
        end
        handles.MultiClass.tbl_cont.array = cell2mat(handles.MultiClass.tbl_cont.array);
        handles.MultiClass.tbl_cont.rownames(remind) = [];
    end
elseif isfield(GDdims,'Regr') % Regression Model 
    handles.labels                                  = Label(I,handles.curlabel);
    handles.modeflag                                = 'regression';
    handles.nclass                                  = 1;
    handles.Regr                                    = GDdims.Regr;
    handles.Regr.labels                             = Label(I,handles.curlabel);
    handles.Regr.cases                              = handles.subjects(I);
    handles.Regr.index_predictions                  = GDdims.Regr.index_predictions(I,handles.curlabel);
    handles.Regr.mean_predictions                   = GDdims.Regr.mean_predictions(I,handles.curlabel);
    handles.Regr.std_predictions                    = GDdims.Regr.std_predictions(I,handles.curlabel );
    handles.Regr.CI1_predictions                    = GDdims.Regr.CI1_predictions(I,handles.curlabel);
    handles.Regr.CI2_predictions                    = GDdims.Regr.CI2_predictions(I,handles.curlabel);
    handles.Regr.best_TR                            = GDdims.bestTR{1}(:,:, handles.curlabel);
    handles.Regr.best_TS                            = GDdims.bestTS{1}(:,:, handles.curlabel);
    handles.Regr.best_P                             = GDdims.bestP{1}(:,handles.curlabel);
    handles.Regr.best_Ppos                          = GDdims.bestPpos{1}(:,handles.curlabel);
    handles.Regr.best_CVperf                        = GDdims.best_CVperf(1);
    handles.Regr.best_TSperf                        = GDdims.best_TSperf(1);
    handles.Regr.best_Complexity                    = GDdims.bestComplexity{1};
    handles.Regr.best_Error                         = GDdims.bestError{1};
    handles.Regr.tbl.colnames                       = {'Cases', 'EXP_LABEL', 'Mean_Pred', 'Std_Pred' };
    handles.Regr.tbl.rownames                       = handles.subjects(handles.Regr.index_predictions);
    handles.Regr.tbl.array                          = [handles.Regr.labels ...
                                                        handles.Regr.mean_predictions ...
                                                        handles.Regr.std_predictions ];
    handles.Regr.tbl_cont.R2                        = handles.Regr.R2(handles.curlabel);
    handles.Regr.tbl_cont.r                         = handles.Regr.r(handles.curlabel);
    handles.Regr.tbl_cont.r_95CI_low                = handles.Regr.r_95CI_low(handles.curlabel);
    handles.Regr.tbl_cont.r_95CI_up                 = handles.Regr.r_95CI_up(handles.curlabel);
    handles.Regr.tbl_cont.t                         = handles.Regr.t(handles.curlabel);
    handles.Regr.tbl_cont.p                         = handles.Regr.p(handles.curlabel);
    handles.Regr.tbl_cont.MAE                       = handles.Regr.MAE(handles.curlabel);
    handles.Regr.tbl_cont.NRSMD                     = handles.Regr.NRSMD(handles.curlabel);
    arr                                             = cell2mat(struct2cell( handles.Regr.tbl_cont));
    handles.Regr.tbl_cont.rownames                  = fieldnames(handles.Regr.tbl_cont);
    handles.Regr.tbl_cont.colnames                  = {'Metric', 'Value'};
    handles.Regr.tbl_cont.array                     = arr;
    
elseif isfield(GDdims,'OneClass') % One-class model
    handles.modeflag                                = 'oneclass';
end