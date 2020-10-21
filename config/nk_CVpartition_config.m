% =========================================================================
% FORMAT res = nk_CVpartition_config(res)
% =========================================================================
%
% Setup repeated nested cross-validation structure 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NeuroMiner 1.0, (c) Nikolaos Koutsouleris 08/2018

function act = nk_CVpartition_config(defaultsfl)

global NM

if ~exist('defaultsfl','var') || isempty(defaultsfl),  defaultsfl = 0; end;

if ~defaultsfl

    CV2frm = 1;
    CV1ps = 'not defined'; CV1fs = CV1ps; CV2ps = CV1fs; CV2fs = CV2ps; 
    CV1pn = 10; CV1fn = 10; CV2pn = 10; CV2fn = 10; 
    Decomp = '';
    CV2ps = ['Define no. of CV2 permutations [ P2 = ' CV2ps ' ]|'];
    CV1ps = ['Define no. of CV1 permutations [ P1 = ' CV1ps ' ]|'];
    
    % Check whether cv structure already exists
    if isfield(NM,'cv')
        if iscell(NM.cv)
            nk_PrintLogo
            ncv = length(NM.cv);
            fprintf('\nMultiple CV structures detected:')
            fprintf('\n================================')
            for i=1:ncv
                % determine size of partitions
                fprintf('\n%g:\t CV2: [%g, &g], CV1: [%g, %g]', ...
                    i, size(NM.cv{i}.TrainInd,1), size(NM.cv{i}.TrainInd,2), ...
                    size(NM.cv{i}.cvin{1,1}.TrainInd,1), size(NM.cv{i}.cvin{1,1}.TrainInd,2))
            end
            actcv = nk_input('What to do',0,'m','Modify existing CV|Remove existing CV|Add new CV|<< Back',1:4);
            switch actcv
                case 1
                    selcv = nk_input(['Select CV for modification [1-' num2str(ncv) ']'],0,e);
                    tcv = NM.cv{selcv};
                case 2
                    selcv = nk_input(['Select CV for removal [1-' num2str(ncv) ']'],0,e);
                    tcv = NM.cv; NM = rmfield(NM,'cv'); cnt=1;
                    for i=1:ncv
                        if i==selcv, continue, end
                        NM.cv{cnt} = tcv{i};
                        cnt=cnt+1;
                    end
                    NM = nk_CVpartition_config(NM);
                case 3
                    selcv = ncv+1; tcv = [];
                case 4
                    return
            end

        else
            ncv = 1;
        end
    end
    fl = true;
    if strcmp(NM.modeflag, 'classification') && isfield(NM,'groupnames') && length(NM.groupnames) > 2
        Decomps = 'not defined';
        if isfield(NM,'TrainParam') && isfield(NM.TrainParam,'RAND') && isfield(NM.TrainParam.RAND,'Decompose')
            switch NM.TrainParam.RAND.Decompose
                case 1
                    Decomps = 'One-vs-One';
                case 2
                    Decomps = 'One-vs-All';
                case 9
                    Decomps = 'Multi-group (no decomposition)';
            end
            Decompn = NM.TrainParam.RAND.Decompose;
        else
            Decompn = 1; fl = false;
        end
        Decomp = ['Define decomposition mode: ' Decomps '|'];
        MenuVec = [11,1:5];
    else
        NM.TrainParam.RAND.Decompose = 1; Decompn = 1;
        MenuVec = [11,1:4];
    end
    
    buildstr = ''; savestr = ''; loadstr = ''; MenuRem = []; CV2prx = '' ; CV1prx = ''; histstr = '';

    %% Define menu options for cross-validation setup
    if isfield(NM,'TrainParam')

        if isfield(NM.TrainParam,'RAND')
            
            if isfield(NM.TrainParam.RAND,'lgoflp'), lgoflp = NM.TrainParam.RAND.lgoflp; end
            CV2reps = ''; MenuVec_extra = [];
            CV2STR_FRAME = '(Pooled) cross-validation';
            if isfield(NM.TrainParam.RAND,'CV2Frame')
               CV2frm = NM.TrainParam.RAND.CV2Frame;
               
               switch NM.TrainParam.RAND.CV2Frame
                    case 2
                       if NM.TrainParam.RAND.OuterPerm>1
                           CV2STR_FRAME = 'Outer Leave-Group-Out/Inner Pooled ';
                       else
                           CV2STR_FRAME = 'Outer Leave-Group-Out/Inner Pooled';
                       end
                       CV2reps = sprintf('Define number of Leave-Group-Out repetitions [ %g repetitions ]|',NM.TrainParam.RAND.OuterPerm);  MenuVec_extra = 12;
                    case 3
                       CV2STR_FRAME = 'Nested Leave-Group-Out';
                       NM.TrainParam.RAND.OuterPerm = 1;
                   case 4
                       CV2STR_FRAME = 'Outer Leave-Group-Out/Inner Leave-Group-In';
                end
            end
            CV2Frame = ['Select cross-validation framework [ ' CV2STR_FRAME ' ]|'];

            if isfield(NM.TrainParam.RAND,'CV2Frame') && NM.TrainParam.RAND.CV2Frame ~= 1
                if ~isfield(NM.TrainParam.RAND,'CV2LCO')
                   CV2STR_LCO = 'not defined';
                   CV2fn = 1;
                else
                   CV2fn = numel(unique(NM.TrainParam.RAND.CV2LCO.ind));
                   CV2STR_LCO = sprintf('%g groups among %g cases',CV2fn, numel(NM.TrainParam.RAND.CV2LCO.ind));
                end
                MenuRem = 3;           
                CV2ps = ['Define Outer (CV2) Leave-Group-Out partitioning [ ' CV2STR_LCO ' ]|'];
                CV2fs = '';
            else
                if isfield(NM.TrainParam.RAND,'OuterPerm')
                    if isfield(NM.TrainParam.RAND,'OuterFold') && ...
                            (NM.TrainParam.RAND.OuterFold == -1 || ...
                            NM.TrainParam.RAND.OuterFold == numel(NM.label))
                        CV2ps = ''; CV2prx = ' [ LOO ]';
                        if strcmp(NM.modeflag,'classification'), NM.SVM.GridParam = 1; end
                        MenuRem = 2;
                    else
                        CV2pn = NM.TrainParam.RAND.OuterPerm; CV2ps = num2str(CV2pn);
                        CV2ps = ['Define no. of Outer (CV2) permutations [ P2 = ' CV2ps ' ]|'];
                    end
                else
                    fl = false;
                end
                if isfield(NM.TrainParam.RAND,'OuterFold')
                    CV2fn = NM.TrainParam.RAND.OuterFold; CV2fs = num2str(CV2fn);
                else
                    fl = false;
                end
                CV2fs = [ 'Define no. of Outer (CV2) folds [ K2 = ' CV2fs CV2prx ' ]|' ];
            end

            if isfield(NM.TrainParam.RAND,'InnerPerm')
                if isfield(NM.TrainParam.RAND,'CV2Frame') && NM.TrainParam.RAND.CV2Frame >= 3 
                    if ~isfield(NM.TrainParam.RAND,'CV1LCO')
                       CV1STR_LCO = 'not defined';
                       CV1fn = 1;
                    else
                       CV1fn = numel(unique(NM.TrainParam.RAND.CV1LCO.ind));
                       CV1STR_LCO = sprintf('%g groups among %g cases',CV1fn, numel(NM.TrainParam.RAND.CV1LCO.ind));
                    end
                    MenuRem = [3 5];    
                    switch NM.TrainParam.RAND.CV2Frame
                        case 3
                            CV1ps = ['Define Inner (CV1) Leave-Group-Out partitioning [ ' CV1STR_LCO ' ]|'];
                        case 4
                            CV1ps = ['Define Inner (CV1) Leave-Group-In partitioning [ ' CV1STR_LCO ' ]|'];
                    end
                    CV1fs = '';
                else
                    if isfield(NM.TrainParam.RAND,'InnerFold') && ... 
                            (NM.TrainParam.RAND.InnerFold == -1 || ...
                            NM.TrainParam.RAND.InnerFold == numel(NM.label) - floor(numel(NM.label) / CV2fn))
                        CV1ps = ''; CV1prx = ' [ LOO ]';
                         if strcmp(NM.modeflag,'classification'), NM.SVM.GridParam = 1; end
                        MenuRem = [MenuRem 4];
                    else
                        CV1pn = NM.TrainParam.RAND.InnerPerm; CV1ps = num2str(CV1pn);
                        CV1ps = ['Define no. of Inner (CV1) permutations [ P1 = ' CV1ps ' ]|'];
                    end
                    CV1fs = ['Define no. of Inner (CV1) folds [ K1 = ' num2str(NM.TrainParam.RAND.InnerFold) CV1prx ' ]|' ];
                end
            else
                fl = false;
            end
            if ~isfield(NM.TrainParam.RAND,'InnerFold') 
                fl = false;
            end
            if fl
                enabled = 2; enabledstr = 'no';
                addremoved2test = 2; addremoved2teststr = ''; 
                if isfield(NM.TrainParam,'RAND') && isfield(NM.TrainParam.RAND,'Eq') 
                    if isfield(NM.TrainParam.RAND.Eq,'enabled') && NM.TrainParam.RAND.Eq.enabled ==1, enabled = 1; enabledstr = 'yes'; end
                    if enabled == 1 && (isfield(NM.TrainParam.RAND.Eq,'addremoved2test') && NM.TrainParam.RAND.Eq.addremoved2test==1), 
                            addremoved2test = 1; addremoved2teststr = ', shuffle to CV1 data'; 
                    end
                end
                switch NM.modeflag
                    case 'regression'
                        histstr = ['Equalize label histogram at the CV1 cycle by undersampling [ ' enabledstr addremoved2teststr ' ]|']; 
                    case 'classification'
                        histstr = ['Equalize class sizes at the CV1 cycle by undersampling [ ' enabledstr addremoved2teststr ' ]|']; 
                end
                MenuVec = [ MenuVec 10];
                buildstr = 'Build Cross-Validation structure|';
                loadstr  = 'Load Cross-Validation structure|';
                MenuVec = [MenuVec 6 7];
                if isfield(NM,'cv')
                    savestr  = 'Save Cross-Validation structure|';
                    MenuVec  = [MenuVec 8];
                end
            end
        end
    end

    if ~isempty(CV2reps), MenuVec = [MenuVec(1:3) MenuVec_extra MenuVec(4:end) ]; end
    if ~isempty(MenuRem), MenuVec(MenuRem) = []; end
    
    nk_PrintLogo

    cprintf('blue','\n****************************************************************************************')
    cprintf('blue','\nSelect appropriate Kx (x=1/2) of folds for your dataset, where:')
    cprintf('blue','\n\tKx < # of available subjects defines K-fold repeated, stratified cross-validation')
    cprintf('blue','\n\tKx = -1 defines LOO cross-validation')
    %cprintf('blue','\n\tKx = -2 defines LOO-per-group cross-validation')
    cprintf('blue','\nIf you choose Kx = -1 no CV permutations will be available.')
    cprintf('blue','\n****************************************************************************************\n')

    act = nk_input('MAIN INTERFACE >> DEFINE PARAMETERS >> CROSS-VALIDATION SETTINGS',0,'mq', ...
                    [CV2Frame ...
                     CV2ps ...
                     CV2fs ...
                     CV2reps ...
                     CV1ps ...
                     CV1fs ...
                     Decomp ...
                     histstr ...
                     buildstr loadstr savestr], MenuVec,1);

    switch act
         case 11
             NM.TrainParam.RAND.CV2Frame = nk_input('Select cross-validation framework',0,'m', ... 
                 'Pooled|Outer Leave-Group-Out/Inner pooled|Nested Leave-Group-Out|Outer Leave-Group-Out/Inner Leave-Group-In',1:4,CV2frm);
             switch NM.TrainParam.RAND.CV2Frame
                 case 1
                    if isfield(NM.TrainParam.RAND,'CV2LCO'), NM.TrainParam.RAND = rmfield(NM.TrainParam.RAND,'CV2LCO'); end
                    if isfield(NM.TrainParam.RAND,'CV1LCO'), NM.TrainParam.RAND = rmfield(NM.TrainParam.RAND,'CV1LCO'); end
                 case 2
                    if isfield(NM.TrainParam.RAND,'CV1LCO'), NM.TrainParam.RAND = rmfield(NM.TrainParam.RAND,'CV1LCO'); end
             end
         case 1
             if isfield(NM.TrainParam.RAND,'CV2Frame') 
                 switch NM.TrainParam.RAND.CV2Frame
                     case 1
                        NM.TrainParam.RAND.OuterPerm = ...
                            nk_input('Number of permutations for Outer (CV2) cross-validation',0,'w1',CV2pn);
                     case {2,3,4}
                         NM.TrainParam.RAND.CV2LCO.ind = ...
                             nk_input('Define index vector for Outer (CV2) case-to-group assignment',0,'i',[],[numel(NM.cases),1]);
                         if NM.TrainParam.RAND.CV2Frame~=2
                             NM.TrainParam.RAND.OuterPerm = 1;
                         end
                         NM.TrainParam.RAND.OuterFold = numel(unique(NM.TrainParam.RAND.CV2LCO.ind));
                 end
             else
                 NM.TrainParam.RAND.OuterPerm = ...
                            nk_input('Number of permutations for Outer (CV2) cross-validation ',0,'w1',CV2pn);
             end
         case 2
             OuterFold = nk_input('Number of folds for Outer (CV2) cross-validation ',0,'i',CV2fn);
             if OuterFold <= 1 && OuterFold ~= -1,
                 errordlg('Specify at least two Outer (CV2) folds for k-fold cross-validation!');
             else
                 NM.TrainParam.RAND.OuterFold = OuterFold;
             end
         case 3
             if isfield(NM.TrainParam.RAND,'CV2Frame') 
                 switch NM.TrainParam.RAND.CV2Frame
                    case {1,2}
                        NM.TrainParam.RAND.InnerPerm = ...
                            nk_input('Number of permutations for Inner (CV1) cross-validation',0,'w1',CV1pn);
                     case {3,4}
                        if NM.TrainParam.RAND.CV2Frame == 3
                              NM.TrainParam.RAND.CV1LCO.ind = ...
                                     nk_input('Define index vector for Inner (CV1) case-to-group assignment',0,'i',[],[numel(NM.cases),1]);
                                 NM.TrainParam.RAND.InnerPerm = 1;
                                 NM.TrainParam.RAND.InnerFold = numel(unique(NM.TrainParam.RAND.CV1LCO.ind));
                        end
                 end
             else
                  NM.TrainParam.RAND.InnerPerm = ...
                            nk_input('Number of permutations for Inner (CV1) cross-validation',0,'w1',CV1pn);
             end
             
         case 4
             InnerFold = nk_input('Number of folds for Inner cross-validation (CV1)',0,'i',CV1fn);
             if InnerFold <= 1 && InnerFold ~= -1,
                 errordlg('Specify at least two CV1 folds for k-fold cross-validation!');
             else
                 NM.TrainParam.RAND.InnerFold = InnerFold;
             end

         case 10

             switch NM.modeflag
                 case 'regression'
                     NM.TrainParam.RAND.Eq.enabled = ...
                         nk_input('Enable histogram equalization at the CV1 level',0,'yes|no',[1,0], enabled);
                 case 'classification'
                     NM.TrainParam.RAND.Eq.enabled = ...
                         nk_input('Enable class size balancing at the CV1 level',0,'yes|no',[1,0], enabled);
             end
             if NM.TrainParam.RAND.Eq.enabled

                 if isfield(NM,'covars') && ~isempty(NM.covars)
                    eqtarget = nk_input('Perform histogram euqalization using the target label or some other covariate', ...
                                        0,'m','Target label|Covariate',[0,1],1);
                    if eqtarget
                         covind = nk_SelectCovariateIndex(NM);
                         NM.TrainParam.RAND.Eq.Covar = NM.covars(:,covind);
                         varstr = NM.covnames{covind};
                    else
                         NM.TrainParam.RAND.Eq.Covar = NM.label;
                         varstr = 'target label';
                    end
                 else
                     NM.TrainParam.RAND.Eq.Covar = NM.label;
                     varstr = 'target label';
                 end
                 switch NM.modeflag
                     case 'classification'
                         NM.TrainParam.RAND.Eq.posnegrat = nk_input('Target ratio bigger / smaller class (1 => classes have sample size after balancing)',0,'e', 1.5);
                     case 'regression'
                         mincount = 1; maxcount = 1; bincount = 10; addremoved2test = 1;
                         if isfield(NM.TrainParam.RAND.Eq,'maxcount'), ...
                             maxcount = NM.TrainParam.RAND.Eq.maxcount; end
                         if isfield(NM.TrainParam.RAND.Eq,'mincount'), ...
                             mincount = NM.TrainParam.RAND.Eq.mincount; end
                         if isfield(NM.TrainParam.RAND.Eq,'bincount'), ...
                             bincount = NM.TrainParam.RAND.Eq.bincount; end
                         fprintf('\n\n************************************************************')
                         fprintf('\nSetup for histogram equalization during CV1 cycle generation')
                         fprintf('\n************************************************************\n')
                         figure;histogram(NM.TrainParam.RAND.Eq.Covar); 
                            ylabel('# of observations in bins'); 
                            xlabel(varstr); 
                            title(['Histogram analysis of ' varstr]);
                         NM.TrainParam.RAND.Eq.maxcount = ...
                             nk_input('Minimum # of observation at the upper end of label histogram',0,'e', maxcount);
                         NM.TrainParam.RAND.Eq.mincount = ...
                             nk_input('Minimum # of observation at the lower end of label histogram',0,'e', mincount);
                         NM.TrainParam.RAND.Eq.bincount = ...
                             nk_input('# Bins for histogram analysis',0,'e',bincount);
                 end
                 NM.TrainParam.RAND.Eq.addremoved2test = ...
                     nk_input('Add removed training subjects to respective CV1 test folds?',0,'yes|no',[1,0],addremoved2test);
             end
         case 5
             if NM.TrainParam.RAND.InnerFold == -1 || NM.TrainParam.RAND.OuterFold == -1
                 NM.TrainParam.RAND.Decompose = ...
                     nk_input('Multi-group decomposition method',0,'m', ...
                     'One-vs-All',2);
             else
                 NM.TrainParam.RAND.Decompose = ...
                     nk_input('Multi-group decomposition method',0,'m', ...
                     'One-vs-One|One-vs-All',[1,2],Decompn);
             end
         case 6
             groups = []; appendfl=false; oldcv=zeros(0,size(NM.label,2));

             % Checkwhether to overwrite or append to current CV structure
             if isfield(NM,'cv')
                 if NM.TrainParam.RAND.OuterFold == size(NM.cv(1).TrainInd,2) && ...
                    NM.TrainParam.RAND.InnerFold == size(NM.cv(1).cvin{1,1}.TrainInd,2)
                        appendfl = nk_input('Append to existing CV structure',0,'m', ...
                            ['No|' ...
                            'Append to Outer CV cycle|' ...
                            'Append to Inner CV cycles|' ...
                            'Rebuild CV structure without new permutations'],0:3,1);
                        if appendfl, oldcv = NM.cv; end
                 end
             end
             if isfield(NM,'groupnames')
                 groupnames = NM.groupnames;
             else
                 groupnames = [];
             end
             if strcmp(NM.modeflag,'classification')
                 for i=1:size(NM.label,2)
                     if ~isempty(oldcv), ioldcv = oldcv(i); else ioldcv=[]; end
                     if size(NM.label,2)>1, grpn = groupnames{i}; else grpn = groupnames; end
                     cv(i) =  nk_MakeCrossFolds(NM.label(:,i), NM.TrainParam.RAND, NM.modeflag, groups, grpn, ioldcv, appendfl);
                 end
             else
                 cv = nk_MakeCrossFolds(NM.label, NM.TrainParam.RAND, NM.modeflag, groups, groupnames, oldcv, appendfl);
             end
             if ~isempty(cv)
                 for i=1:numel(cv)
                     switch appendfl
                         case 1
                            NM.cv(i).TrainInd = [NM.cv(i).TrainInd; cv(i).TrainInd];
                            NM.cv(i).TestInd 	= [NM.cv(i).TestInd; cv(i).TestInd];
                            NM.cv(i).cvin	= [NM.cv(i).cvin; cv(i).cvin];
                            if strcmp(NM.modeflag,'classification')
                                NM.cv(i).class = [NM.cv(i).class; cv(i).class];
                                NM.cv(i).classnew = [NM.cv(i).classnew; cv(i).classnew];
                            end
                         otherwise
                            NM.cv(i) = cv(i);
                     end
                 end
             end
         case 7
             loadcv = nk_FileSelector(1,'matrix','Select cross-validation structure', 'CVstruct_.*\.mat');
             load(loadcv,'cv');
             NM.cv = cv;
         case 8
             savecv = nk_input('Path & filename for cross-validation structure',0,'s');
             cv =  NM.cv;
             save(['CVstruct_' savecv '.mat'],'cv');
             
        case 12
             NM.TrainParam.RAND.OuterPerm = nk_input('Number of repetitions for Outer (CV2) cross-validation',0,'w1',CV2pn);
                             
    end
else
    NM.TrainParam.RAND.CV2Frame     = 1;
    NM.TrainParam.RAND.OuterFold    = 10;
    NM.TrainParam.RAND.OuterPerm    = 1;
    NM.TrainParam.RAND.InnerFold    = 10;
    NM.TrainParam.RAND.InnerPerm    = 1;
end

