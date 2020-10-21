function [InputParam, TrainedParam, SrcParam] = nk_GenPreprocSequence(InputParam, TemplParam, SrcParam, TrainedParam)

global MODEFL NM VERBOSE

if isfield(TemplParam,'ACTPARAM')
    
    % Prepare preprocessing
    lact                    = numel(TemplParam.ACTPARAM);
    InputParam.P            = cell(lact,1);
    InputParam.CV1perm      = SrcParam.CV1perm;
    InputParam.CV1fold      = SrcParam.CV1fold;
    actionseq               = cell(1,lact);
    InputParam.curclass     = SrcParam.u;
  
    % Loop through ACTPARAM sequence
    for ac = 1:lact
        
        actionseq{ac} = TemplParam.ACTPARAM{ac}.cmd;
        
        switch actionseq{ac}
            
            case 'labelimpute'
                
                if VERBOSE, fprintf('\n* Impute NaNs in label vector.'); end
                InputParam.P{ac}.LABELIMPUTE = TemplParam.ACTPARAM{ac}.LABELIMPUTE;
                SrcParam.NaNflag = true;
                InputParam.P{ac}.BINMOD = TemplParam.BINMOD;
                    
            case 'impute'
                
                if VERBOSE, fprintf('\n* Impute NaNs in matrix block of %g columns.',sum(TemplParam.ACTPARAM{ac}.IMPUTE.blockind)); end
                InputParam.P{ac}.IMPUTE = TemplParam.ACTPARAM{ac}.IMPUTE;
        
            case 'copyfeat'

                if VERBOSE, fprintf('\n* Use all data.'); end
            
            case 'normalize'
                
                if ~isempty(SrcParam.covars)
                    InputParam.P{ac}.IND = TemplParam.ACTPARAM{ac}.IND;
                    InputParam.P{ac}.TsInd = [];
                    if VERBOSE, 
                        fprintf('\n* PER-GROUP NORMALIZATION AND ZERO-VARIANCE REMOVAL ACROSS GROUPS')
                        fprintf('\n\t- Normalize data to mean effects in group: %s', NM.covnames{TemplParam.ACTPARAM{ac}.IND});
                    end
                    if isfield(SrcParam,'TrX'),         InputParam.P{ac}.TrInd        = SrcParam.covars( SrcParam.TrX, TemplParam.ACTPARAM{ac}.IND );   end
                    if isfield(SrcParam,'TrI'),         InputParam.P{ac}.TsInd{end+1} = SrcParam.covars( SrcParam.TrI, TemplParam.ACTPARAM{ac}.IND );   end
                    if isfield(SrcParam,'CVI'),         InputParam.P{ac}.TsInd{end+1} = SrcParam.covars( SrcParam.CVI, TemplParam.ACTPARAM{ac}.IND );   end
                    if isfield(SrcParam,'TsI'),         InputParam.P{ac}.TsInd{end+1} = SrcParam.covars( SrcParam.TsI, TemplParam.ACTPARAM{ac}.IND );   end
                    if ~isempty(SrcParam.covars_oocv),  InputParam.P{ac}.TsInd{end+1} = SrcParam.covars_oocv( :, TemplParam.ACTPARAM{ac}.IND );   end
                    
                    if ~isempty(SrcParam.iTr), 
                        InputParam.P{ac}.TrInd(SrcParam.iTr,:)=[]; 
                        InputParam.P{ac}.TsInd{1}(SrcParam.iTr,:)=[];
                    end
                    if ~isempty(SrcParam.iCV), InputParam.P{ac}.TsInd{2}(SrcParam.iCV,:)=[]; end
                    if ~isempty(SrcParam.iTs), InputParam.P{ac}.TsInd{3}(SrcParam.iTs,:)=[]; end
                end
                
            case 'correctnuis'
                
                if ~isempty(SrcParam.covars)
                    InputParam.P{ac}.COVAR = TemplParam.ACTPARAM{ac}.COVAR;
                    if VERBOSE, 
                        fprintf('\n* ADJUSTING DATA FOR COVARIATE EFFECTS')
                        fprintf('\n\t- Nuisance covariate: %s', NM.covnames{TemplParam.ACTPARAM{ac}.COVAR})
                    end
                    InputParam.P{ac}.TsCovars = [];
                    if isfield(SrcParam,'TrX'),         InputParam.P{ac}.TrCovars        = SrcParam.covars( SrcParam.TrX, TemplParam.ACTPARAM{ac}.COVAR );   end
                    if isfield(SrcParam,'TrI'),         InputParam.P{ac}.TsCovars{end+1} = SrcParam.covars( SrcParam.TrI, TemplParam.ACTPARAM{ac}.COVAR );   end
                    if isfield(SrcParam,'CVI'),         InputParam.P{ac}.TsCovars{end+1} = SrcParam.covars( SrcParam.CVI, TemplParam.ACTPARAM{ac}.COVAR );   end
                    if isfield(SrcParam,'TsI'),         InputParam.P{ac}.TsCovars{end+1} = SrcParam.covars( SrcParam.TsI, TemplParam.ACTPARAM{ac}.COVAR );   end
                    if ~isempty(SrcParam.covars_oocv),  InputParam.P{ac}.TsCovars{end+1} = SrcParam.covars_oocv( :,       TemplParam.ACTPARAM{ac}.COVAR );   end
                    if isfield(TemplParam.ACTPARAM{ac},'METHOD') && TemplParam.ACTPARAM{ac}.METHOD==2
                        if VERBOSE, fprintf('\n\t- Method: Combat'); end
                        InputParam.P{ac}.METHOD = 2; InputParam.P{ac}.TsMod = [];
                        if TemplParam.ACTPARAM{ac}.MCOVARLABEL
                            uL = unique(NM.label); uL(isnan(uL))=[]; nuL = numel(uL);
                            switch NM.modeflag
                                case 'classification'
                                    dummy_labels = nk_MakeDummyVariables(NM.label); 
                                    if nuL==2, dummy_labels(:,2)=[]; end
                                case 'regression'
                                    dummy_labels = NM.label;
                            end
                            covars = [ SrcParam.covars(:, TemplParam.ACTPARAM{ac}.MCOVAR) dummy_labels ];
                            if ~isempty(SrcParam.covars_oocv) 
                                covars_oocv = [ SrcParam.covars_oocv(:,TemplParam.ACTPARAM{ac}.MCOVAR) zeros(size(SrcParam.covars_oocv,1), size(dummy_labels,2)) ];
                            end
                        else
                            covars = SrcParam.covars( : , TemplParam.ACTPARAM{ac}.MCOVAR ); if ~isempty(SrcParam.covars_oocv), covars_oocv = SrcParam.covars_oocv ( : , TemplParam.ACTPARAM{ac}.MCOVAR ); end
                        end
                        if isfield(SrcParam,'TrX'),         InputParam.P{ac}.TrMod        = covars( SrcParam.TrX, :);   end
                        if isfield(SrcParam,'TrI'),         InputParam.P{ac}.TsMod{end+1} = covars( SrcParam.TrI, :);   end
                        if isfield(SrcParam,'CVI'),         InputParam.P{ac}.TsMod{end+1} = covars( SrcParam.CVI, :);   end
                        if isfield(SrcParam,'TsI'),         InputParam.P{ac}.TsMod{end+1} = covars( SrcParam.TsI, :);   end
                        if ~isempty(SrcParam.covars_oocv),  InputParam.P{ac}.TsMod{end+1} = covars_oocv;   end
                    else
                        if VERBOSE, fprintf('\n\t- Method: Partial correlations analysis'); end
                        InputParam.P{ac}.METHOD = 1;
                    end
                    if ~isempty(SrcParam.iTr), 
                        InputParam.P{ac}.TrCovars(SrcParam.iTrX,:)   = []; 
                        InputParam.P{ac}.TsCovars{1}(SrcParam.iTr,:) = []; 
                        if InputParam.P{ac}.METHOD == 2
                            InputParam.P{ac}.TrMod(SrcParam.iTrX,:)   = []; 
                            InputParam.P{ac}.TsMod{1}(SrcParam.iTr,:) = []; 
                        end
                    end
                    if ~isempty(SrcParam.iCV), InputParam.P{ac}.TsCovars{2}(SrcParam.iCV,:)=[]; end
                    if ~isempty(SrcParam.iTs), InputParam.P{ac}.TsCovars{3}(SrcParam.iTs,:)=[]; end
                    if InputParam.P{ac}.METHOD == 2
                        if ~isempty(SrcParam.iCV), InputParam.P{ac}.TsMod{2}(SrcParam.iCV,:)=[]; end
                        if ~isempty(SrcParam.iTs), InputParam.P{ac}.TsMod{3}(SrcParam.iTs,:)=[]; end
                    end
                end
                if InputParam.P{ac}.METHOD == 1;
                    if isfield(TemplParam.ACTPARAM{ac},'INTERCEPT')
                        InputParam.P{ac}.INTERCEPT = TemplParam.ACTPARAM{ac}.INTERCEPT-1;
                        if VERBOSE, 
                            switch InputParam.P{ac}.INTERCEPT
                                case 0
                                    fprintf('\n\t-> Not including intercept.'); 
                                case 1
                                    fprintf('\n\t-> Including intercept.'); 
                            end  
                        end
                    end
                    if isfield(TemplParam.ACTPARAM{ac},'COVDIR')
                        InputParam.P{ac}.COVDIR = TemplParam.ACTPARAM{ac}.COVDIR-1;
                        if VERBOSE, 
                            switch InputParam.P{ac}.COVDIR
                                case 0
                                    fprintf('\n\t-> Covariate effects will be removed from data.'); 
                                case 1
                                    fprintf('\n\t-> Covariate effects will be increased in data.')
                            end  
                        end
                    end
                    if isfield(TemplParam.ACTPARAM{ac},'BETAEXT') && ~isempty(TemplParam.ACTPARAM{ac}.BETAEXT)
                        InputParam.P{ac}.BETAEXT = TemplParam.ACTPARAM{ac}.BETAEXT;
                        if VERBOSE,fprintf('\n\t Beta parameter(s) computed in an OOT-sample will be used.'); end
                    else
                        if isfield(TemplParam.ACTPARAM{ac},'SUBGROUP') && ~isempty(TemplParam.ACTPARAM{ac}.SUBGROUP)
                            InputParam.P{ac}.SUBGROUP = TemplParam.ACTPARAM{ac}.SUBGROUP(SrcParam.TrX,:);
                            if VERBOSE,fprintf('\n\t-> Beta parameter(s) will be computed from a specific subgroup.'); end
                        end
                    end
                else
                    if isfield(TemplParam.ACTPARAM{ac},'SUBGROUP') && ~isempty(TemplParam.ACTPARAM{ac}.SUBGROUP)
                        InputParam.P{ac}.SUBGROUP = TemplParam.ACTPARAM{ac}.SUBGROUP(SrcParam.TrX,:);
                        if VERBOSE,fprintf('\n\t-> Combat parameter(s) will be computed from a specific subgroup.'); end
                    end
                     InputParam.P{ac}.COVDIR=0;
                     InputParam.P{ac}.INTERCEPT=0;
                end
                
            case 'remmeandiff'
                
                if isfield(NM,'covars') && ~isempty(NM.covars)
                    InputParam.P{ac}.sIND = TemplParam.ACTPARAM{ac}.sIND; InputParam.P{ac}.dIND = TemplParam.ACTPARAM{ac}.dIND;
                    if VERBOSE,
                        fprintf('\n* OFFSET CORRECTION')
                        fprintf('\n\t - Compute group mean offset of : %s', NM.covnames{TemplParam.ACTPARAM{ac}.sIND})
                        fprintf('\n\t - Remove global mean and offsets from : %s', NM.covnames{TemplParam.ACTPARAM{ac}.dIND})
                    end
                    if isfield(SrcParam,'TrX'),         
                        InputParam.P{ac}.sTrInd        = SrcParam.covars( SrcParam.TrX, TemplParam.ACTPARAM{ac}.sIND );
                        InputParam.P{ac}.dTrInd        = SrcParam.covars( SrcParam.TrX, TemplParam.ACTPARAM{ac}.dIND );
                        if ~isempty(SrcParam.iTrX), 
                            InputParam.P{ac}.sTrInd( SrcParam.iTrX,:) = []; 
                            InputParam.P{ac}.dTrInd( SrcParam.iTrX,:) = []; 
                        end
                    end
                    if isfield(SrcParam,'TrI'),         
                        InputParam.P{ac}.sTsInd{1}     = SrcParam.covars( SrcParam.TrI, TemplParam.ACTPARAM{ac}.sIND );   
                        InputParam.P{ac}.dTsInd{1}     = SrcParam.covars( SrcParam.TrI, TemplParam.ACTPARAM{ac}.dIND );
                        if ~isempty(SrcParam.iTr), 
                            InputParam.P{ac}.sTsInd{1}( SrcParam.iTr,:) = []; 
                            InputParam.P{ac}.dTsInd{1}( SrcParam.iTr,:) = []; 
                        end
                    end
                    if isfield(SrcParam,'CVI'),         
                        InputParam.P{ac}.sTsInd{end+1} = SrcParam.covars( SrcParam.CVI, TemplParam.ACTPARAM{ac}.sIND );
                        InputParam.P{ac}.dTsInd{end+1} = SrcParam.covars( SrcParam.CVI, TemplParam.ACTPARAM{ac}.dIND );
                        if ~isempty(SrcParam.iCV), 
                            InputParam.P{ac}.sTsInd{2}( SrcParam.iCV,:) = []; 
                            InputParam.P{ac}.dTsInd{2}( SrcParam.iCV,:) = []; 
                        end
                    end
                    if isfield(SrcParam,'TsI'),         
                        InputParam.P{ac}.sTsInd{end+1} = SrcParam.covars( SrcParam.TsI, TemplParam.ACTPARAM{ac}.sIND );
                        InputParam.P{ac}.dTsInd{end+1} = SrcParam.covars( SrcParam.TsI, TemplParam.ACTPARAM{ac}.dIND );
                        if ~isempty(SrcParam.iTs), 
                            InputParam.P{ac}.sTsInd{3}( SrcParam.iTs,:) = []; 
                            InputParam.P{ac}.dTsInd{3}( SrcParam.iTs,:) = []; 
                        end
                    end
                    if ~isempty(SrcParam.covars_oocv)
                        InputParam.P{ac}.sTsInd{end+1} = SrcParam.covars_oocv( : ,    TemplParam.ACTPARAM{ac}.sIND );   
                        InputParam.P{ac}.dTsInd{end+1} = SrcParam.covars_oocv( : ,    TemplParam.ACTPARAM{ac}.dIND );
                    end
                    
                end
                
            case 'standardize' 
                
                % **************** STANDARDIZATION ***************
               
                if VERBOSE, fprintf('\n* Z-NORMALIZATION'); end
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX.opt)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                end
                if isfield(TemplParam.ACTPARAM{ac},'METHOD') && ~isempty(TemplParam.ACTPARAM{ac}.METHOD)
                    InputParam.P{ac}.method = TemplParam.ACTPARAM{ac}.METHOD;
                end
                if isfield(TemplParam.ACTPARAM{ac},'sIND') && ~isempty(TemplParam.ACTPARAM{ac}.sIND)
                    InputParam.P{ac}.sIND = TemplParam.ACTPARAM{ac}.sIND ;
                    if isfield(TemplParam.ACTPARAM{ac},'CALIBUSE') && TemplParam.ACTPARAM{ac}.CALIBUSE
                        InputParam.P{ac}.CALIBUSE      = true;
                    else
                        InputParam.P{ac}.CALIBUSE      = false;
                    end
                    if isfield(SrcParam,'TrX'),  
                        InputParam.P{ac}.sTrInd        = SrcParam.covars( SrcParam.TrX, TemplParam.ACTPARAM{ac}.sIND );
                        if ~isempty(SrcParam.iTrX), InputParam.P{ac}.sTrInd(SrcParam.iTrX,:)=[]; end
                        denom = numel(InputParam.P{ac}.sTrInd);
                    elseif isfield(SrcParam,'TrI')
                        InputParam.P{ac}.sTsInd{1}     = SrcParam.covars( SrcParam.TrI, TemplParam.ACTPARAM{ac}.sIND );
                        if ~isempty(SrcParam.iTr), InputParam.P{ac}.sTsInd{1}(SrcParam.iTr,:)=[]; end
                        denom = numel(InputParam.P{ac}.sTsInd{1});
                    end
                     if VERBOSE,fprintf('\n\tMean & SD computation restricted to N=%g (%1.0f%% of training sample).', ...
                        sum(InputParam.P{ac}.sTrInd), sum(InputParam.P{ac}.sTrInd) * 100 / denom); end
                else
                    InputParam.P{ac}.sTrInd = [];
                end
                
                %InputParam.P{ac}.dTsInd = cell(4,1);
                if isfield(TemplParam.ACTPARAM{ac},'dIND') && ~isempty(TemplParam.ACTPARAM{ac}.dIND)  
                    InputParam.P{ac}.dIND = TemplParam.ACTPARAM{ac}.dIND ;
                    if isfield(SrcParam,'TrX'),         
                        InputParam.P{ac}.dTrInd    = SrcParam.covars( SrcParam.TrX, TemplParam.ACTPARAM{ac}.dIND );
                        if ~isempty(SrcParam.iTrX), InputParam.P{ac}.dTrInd(SrcParam.iTrX,:)=[]; end
                    end
                    if isfield(SrcParam,'TrI'),         
                        InputParam.P{ac}.dTsInd{1} = SrcParam.covars( SrcParam.TrI, TemplParam.ACTPARAM{ac}.dIND );
                        if ~isempty(SrcParam.iTr), InputParam.P{ac}.dTsInd{1}(SrcParam.iTr,:)=[]; end
                    end
                    if isfield(SrcParam,'CVI'),         
                        InputParam.P{ac}.dTsInd{2} = SrcParam.covars( SrcParam.CVI, TemplParam.ACTPARAM{ac}.dIND );
                        if ~isempty(SrcParam.iCV), InputParam.P{ac}.dTsInd{2}(SrcParam.iCV,:)=[]; end
                    end
                    if isfield(SrcParam,'TsI'),         
                        InputParam.P{ac}.dTsInd{3} = SrcParam.covars( SrcParam.TsI, TemplParam.ACTPARAM{ac}.dIND );
                        if ~isempty(SrcParam.iTs), InputParam.P{ac}.dTsInd{3}(SrcParam.iTs,:)=[]; end
                    end
                    if ~isempty(SrcParam.covars_oocv)   
                        InputParam.P{ac}.dTsInd{4} = SrcParam.covars_oocv( : , TemplParam.ACTPARAM{ac}.dIND);
                    end
                end
                
                InputParam.P{ac}.WINSOPT = TemplParam.ACTPARAM{ac}.WINSOPT;
                
            case 'spatialfilter'
                
                if VERBOSE, fprintf('\n* SPATIAL FILTERING'); end
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX.opt)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                end
                InputParam.P{ac}.SPATIAL = TemplParam.ACTPARAM{ac}.SPATIAL;
                
            case 'unitnormalize'
                % **************** UNIT VECTOR NORMALIZATION ***************
                if VERBOSE, fprintf('\n* UNIT VECTOR NORMALIZATION'); end
                InputParam.P{ac}.METHOD = TemplParam.ACTPARAM{ac}.METHOD;
                
            case 'discretize'
                
                % **************** DISCRETIZATION *****************
                if VERBOSE,fprintf('\n* DISCRETIZATION'); end                   
                InputParam.P{ac}.DISCRET = TemplParam.ACTPARAM{ac}.DISCRET;
                InputParam.P{ac}.method = str2func('PerfDiscretizeObj');
                
            case 'symbolize'
                
                % **************** SYMBOLIZATION *****************
                if VERBOSE,fprintf('\n* SYMBOLIZATION'); end               
                InputParam.P{ac}.SYMBOL = TemplParam.ACTPARAM{ac}.SYMBOL;
                InputParam.P{ac}.method = str2func('PerfSymbolizeObj');
                                
            case 'reducedim'
                
                % **************** DIM. REDUCTION *****************
                if ~strcmp(TemplParam.ACTPARAM{ac}.DR.RedMode,'none')
                    if VERBOSE,fprintf('\n* DIMENSIONALITY REDUCTION'); end
                    InputParam.P{ac}.DR = TemplParam.ACTPARAM{ac}.DR;
                    InputParam.P{ac}.DR.Modus = 'JDQR';
                    switch MODEFL
                        case 'classification' 
                            if TemplParam.BINMOD
                                InputParam.P{ac}.DR.labels = SrcParam.BinaryTrainLabel;
                            else
                                InputParam.P{ac}.DR.labels = SrcParam.MultiTrainLabel;
                            end
                        case 'regression'
                            InputParam.P{ac}.DR.labels = SrcParam.TrainLabel;
                    end
                    if strcmp(TemplParam.ACTPARAM{ac}.DR.RedMode,'PLS')
                        InputParam.P{ac}.DR.PLS.V = TemplParam.ACTPARAM{ac}.DR.PLS.V(SrcParam.TrX,:);
                        if sum(SrcParam.iTr), InputParam.P{ac}.DR.PLS.V(SrcParam.iTr,:)=[]; end
                    end
                    
                else
                    InputParam.P{ac}.DR.RedMode = 'none';
                    InputParam.P{ac}.DR.dims = size(InputParam.Tr,2);
                end
                
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                    PX = nk_ReturnParamChain(TemplParam.ACTPARAM{ac});
                    InputParam.P{ac}.DR.Params = PX.Params;
                    InputParam.P{ac}.DR.Params_desc = PX.Params_desc;
                end
                
            case 'scale'
                
                % ******************* SCALING *********************
                if isfield(TemplParam.ACTPARAM{ac}.SCALE,'ZeroOne')
                    if VERBOSE,
                        switch TemplParam.ACTPARAM{ac}.SCALE.ZeroOne
                            case 1
                                fprintf('\n* SCALING [0 1]')
                            case 2
                                fprintf('\n* SCALING [-1 1]')
                        end
                    end
                end
                InputParam.P{ac}.SCALE = TemplParam.ACTPARAM{ac}.SCALE;
           
            case 'rankfeat'
                
                % **************** FEATURE RANKING ***************
                if VERBOSE, fprintf('\n* FEATURE WEIGHTING'); end
                InputParam.P{ac}.RANK = TemplParam.ACTPARAM{ac}.RANK;
                InputParam.P{ac}.RANK.curlabel = TemplParam.ACTPARAM{ac}.RANK.label(SrcParam.TrX,:);
                if ~isempty(SrcParam.iTrX), InputParam.P{ac}.RANK.curlabel(SrcParam.iTrX,:)=[]; end
                if isfield( TemplParam.ACTPARAM{ac}.RANK,'glabel' )
                    % glabel is a logical vector
                    InputParam.P{ac}.RANK.curglabel = TemplParam.ACTPARAM{ac}.RANK.glabel(SrcParam.TrX);
                    if ~isempty(SrcParam.iTrX), InputParam.P{ac}.RANK.curglabel(SrcParam.iTrX,:)=[]; end
                end
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                    PX = nk_ReturnParamChain(TemplParam.ACTPARAM{ac});
                    InputParam.P{ac}.RANK.Params = PX.Params;
                    InputParam.P{ac}.RANK.Params_desc = PX.Params_desc;
                end
                
            case 'extfeat'
                % **************** FEATURE EXTRACTION ***************
                if VERBOSE, fprintf('\n* WEIGHTING-BASED FEATURE GENERATION'); end
                InputParam.P{ac}.W_ACT = TemplParam.ACTPARAM{ac}.W_ACT;
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                    PX = nk_ReturnParamChain(TemplParam.ACTPARAM{ac});
                    InputParam.P{ac}.W_ACT.Params = PX.Params;
                    InputParam.P{ac}.W_ACT.Params_desc = PX.Params_desc;
                end
                
            case 'extdim'
                 % ******************* COMPONENT EXTRACTION (HELPER OF DIM. REDUCTION) *********************
                if VERBOSE, fprintf('\n* COMPONENT EXTRACTION FOLLOWING DIMENSIONALITY REDUCTION'); end
                InputParam.P{ac}.EXTDIM = TemplParam.ACTPARAM{ac}.EXTDIM;
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                    PX = nk_ReturnParamChain(TemplParam.ACTPARAM{ac});
                    InputParam.P{ac}.EXTDIM.Params = PX.Params;
                    InputParam.P{ac}.EXTDIM.Params_desc = PX.Params_desc;
                end
                
            case 'elimzero'
                 % ******************* ELIMINATE ZERO VAR, NAN/INF FEATURES *********************
                if VERBOSE, fprintf('\n* PRUNE ATTRIBUTES'); end
                InputParam.P{ac} =  TemplParam.ACTPARAM{ac};
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX.opt)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                end
                
            case 'remvarcomp'
                 % ******************* REMOVE VARIANCE COMPONENTS  *********************
                if VERBOSE, fprintf('\n* EXTRACT VARIANCE COMPONENTS'); end
                InputParam.P{ac}.REMVARCOMP = TemplParam.ACTPARAM{ac}.REMVARCOMP;
                InputParam.P{ac}.REMVARCOMP.G = TemplParam.ACTPARAM{ac}.REMVARCOMP.G(SrcParam.TrX,:);
                if ~isempty(SrcParam.iTrX), 
                    InputParam.P{ac}.REMVARCOMP.G(SrcParam.iTrX,:)=[]; 
                end
                if isfield(TemplParam.ACTPARAM{ac}.REMVARCOMP,'SUBGROUP') 
                    switch TemplParam.ACTPARAM{ac}.REMVARCOMP.SUBGROUP.flag
                        case 1
                            InputParam.P{ac}.REMVARCOMP.indX = true(numel(find(SrcParam.TrX)),1);
                        case 2
                            InputParam.P{ac}.REMVARCOMP.indX = TemplParam.ACTPARAM{ac}.REMVARCOMP.SUBGROUP.ind(SrcParam.TrX);
                        case 3
                            n = numel(find(SrcParam.TrX))/100*TemplParam.ACTPARAM{ac}.REMVARCOMP.SUBGROUP.indperc;
                            ind = true(numel(find(SrcParam.TrX)),1); randind = randperm(numel(ind),n);
                            InputParam.P{ac}.REMVARCOMP.indX = ind(randind);
                        case 4
                            InputParam.P{ac}.REMVARCOMP.indX = true(numel(InputParam.C,1));
                        case 5
                            InputParam.P{ac}.REMVARCOMP.indX = TemplParam.ACTPARAM{ac}.REMVARCOMP.SUBGROUP.ind;
                        case 6
                            n = size(InputParam.C,1)/100*TemplParam.ACTPARAM{ac}.REMVARCOMP.SUBGROUP.indperc;
                            ind = true(size(InputParam.C,1),1); randind = randperm(numel(ind),n);
                            InputParam.P{ac}.REMVARCOMP.indX = ind(randind);
                    end
                    if ~isempty(SrcParam.iTrX), InputParam.P{ac}.REMVARCOMP.indX(SrcParam.iTrX)=[]; end
                end
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                    PX = nk_ReturnParamChain(TemplParam.ACTPARAM{ac});
                    InputParam.P{ac}.REMVARCOMP.Params = PX.Params;
                    InputParam.P{ac}.REMVARCOMP.Params_desc = PX.Params_desc;
                end
                
            case 'devmap'
                 % ******************* SENSITIZE FEATURES TO DEVIATION FROM NORMATIVE MODEL *********************
                InputParam.P{ac}.DEVMAP = TemplParam.ACTPARAM{ac}.DEVMAP;
              
                if isfield(SrcParam,'TrX'),         
                    InputParam.P{ac}.TrInd    = SrcParam.TrX;
                    if ~isempty(SrcParam.iTrX), InputParam.P{ac}.TrInd(SrcParam.iTrX,:)=[]; end
                end
                if isfield(SrcParam,'TrI'),     
                    InputParam.P{ac}.TsInd{1} = SrcParam.TrI ;
                    if ~isempty(SrcParam.iTr), InputParam.P{ac}.TsInd{1}(SrcParam.iTr,:)=[]; end
                end
                if isfield(SrcParam,'CVI'),         
                    InputParam.P{ac}.TsInd{2} = SrcParam.CVI;
                    if ~isempty(SrcParam.iCV), InputParam.P{ac}.TsInd{2}(SrcParam.iCV,:)=[]; end
                end
                if isfield(SrcParam,'TsI'),         
                    InputParam.P{ac}.TsInd{3} = SrcParam.TsI;
                    if ~isempty(SrcParam.iTs), InputParam.P{ac}.TsInd{3}(SrcParam.iTs,:)=[]; end
                end
                if numel(InputParam.Ts)==4
                    InputParam.P{ac}.TsInd{4} = true(size(InputParam.Ts{4},1),1);
                end
                
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                    PX = nk_ReturnParamChain(TemplParam.ACTPARAM{ac});
                    InputParam.P{ac}.DEVMAP.Params = PX.Params;
                    InputParam.P{ac}.DEVMAP.Params_desc = PX.Params_desc;
                end
        end
    end
    if VERBOSE, fprintf('\nPreprocessing sequence setup completed. Executing ...'); end
    if exist('TrainedParam','var') && ~isempty(TrainedParam)
        [SrcParam, InputParam, TrainedParam] = nk_PerfPreprocessObj_core(SrcParam, InputParam, TrainedParam, actionseq);
    else
        [SrcParam, InputParam, TrainedParam] = nk_PerfPreprocessObj_core(SrcParam, InputParam, [], actionseq);
    end
    clear actionseq
else
    if VERBOSE,fprintf('\nNo preprocessing sequence detected. Aborting ...'); end
    TrainedParam = [];
end

end