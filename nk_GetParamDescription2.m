function vargout = nk_GetParamDescription2(res, params, action, descrip, varind)

global NMinfo

if exist('descrip','var') && ~isempty(descrip)
    vargout = descrip;
else
    vargout = [];
end

switch action
    
    case 'VarDesc'
        
        if isfield(params,'datadescriptor')
            for i=1:numel(varind)
                if isfield(params.datadescriptor{varind(i)},'type')
                    switch params.datadescriptor{varind(i)}.type
                        case 0
                            vargout.datadescriptor{i} = 'MATLAB matrix';
                        case 1
                            vargout.datadescriptor{i} = 'ANALYZE / NIFTI data';
                        case 2
                            vargout.datadescriptor{i} = 'Freesurfer data (MGH / MGZ)';
                        case 3
                            vargout.datadescriptor{i} = '';
                        case 4
                            vargout.datadescriptor{i} = '';
                        case 5
                            vargout.datadescriptor{i} = '';
                    end
                else
                    vargout.datadescriptor{i}= '';
                end
                if isfield(params.datadescriptor{varind(i)},'desc') && ~isempty(params.datadescriptor{varind(i)}.desc)
                    vargout.datadescriptor{i} = [vargout.datadescriptor{i} ' [' params.datadescriptor{varind(i)}.desc ']'];
                end
            end
        else
            vargout.datadescriptor{i} = 'NA';
        end
    
    case 'Version'
        vargout.version = sprintf('%s, %s', NMinfo.info.name, NMinfo.info.ver);
        vargout.author  = sprintf('%s (%s)', NMinfo.info.author, NMinfo.info.email);
        vargout.date    = sprintf('%s', NMinfo.info.datever);        

    % %%%%%%%%%% Preprocessing Parameters %%%%%%%%%%%    
    case 'PreProc'
        
        preprocact = [];
        if ~isempty(params) && isfield(params,'BINMOD')
            switch params.BINMOD
                case 1
                    binmodestr = 'binary';
                case 0
                    binmodestr = 'multi-group';
            end
            groupmode = ['Group processing mode: ' binmodestr];
        else
            groupmode = 'Group processing mode: undefined';
        end
        targetscalestr = 'NA'; labelimputestr = 'NA';
        if isfield(params,'LABELMOD') 
            if isfield(params.LABELMOD,'TARGETSCALE') && ~isempty(params.LABELMOD.TARGETSCALE)
                if params.LABELMOD.TARGETSCALE, 
                    targetscalestr = 'Target scaling [0, 1]';
                else
                    targetscalestr = '';
                end
            end
            if isfield(params.LABELMOD,'POLYNOM') && ~isempty(params.LABELMOD.POLYNOM)
                if numel(targetscalestr)>1, targetscalestr = [targetscalestr, ', ']; end
                targetscalestr = [ targetscalestr 'Exponential transformation: label.^' num2str(params.LABELMOD.POLYNOM) ];
            end
            if isfield(params.LABELMOD,'LABELIMPUTE') && ~isempty(params.LABELMOD.LABELIMPUTE)
                switch params.LABELMOD.LABELIMPUTE.method
                    case 'ml'
                        impstr = 'ML algorithm used for training the models';
                    case 'singlemean'
                        impstr = 'single-subject median replacement';
                    case 'mean'
                        impstr = 'Feature-wise mean replacement';
                    case {'manhattan','dist'}
                        impstr = sprintf('median of %g nearest neighbors (%s)',params.LABELMOD.LABELIMPUTE.k, 'Manhattan');
                    case {'euclidean','dist2'}
                        impstr = sprintf('median of %g nearest neighbors (%s)',params.LABELMOD.LABELIMPUTE.k, 'Euclidean');
                    case 'seuclidean'
                        impstr = sprintf('median of %g nearest neighbors (%s)',params.LABELMOD.LABELIMPUTE.k, 'Seuclidean');
                    case 'cosine'
                        impstr = sprintf('median of %g nearest neighbors (%s)',params.LABELMOD.LABELIMPUTE.k, 'Cosine similarity');
                    case 'mahalanobis'
                        impstr = sprintf('median of %g nearest neighbors (%s)',params.LABELMOD.LABELIMPUTE.k, 'Mahalanobis');
                    case 'jaccard'
                        impstr = sprintf('median of %g nearest neighbors (%s)',params.LABELMOD.LABELIMPUTE.k, 'Jaccard');
                    case 'hamming'
                        impstr = sprintf('median of %g nearest neighbors (%s)',params.LABELMOD.LABELIMPUTE.k, 'Hamming');
                    case 'none'
                        impstr = 'none';
                end
                labelimputestr = sprintf('Label Imputation: %s', impstr);
            end
        end
        spatialfilterstr = 'NA';
        if isfield(params,'SPATIAL') && ~isempty(params.SPATIAL)
            switch params.SPATIAL.cubetype
                case 1
                    spatialfilterstr = 'No spatial filtering';
                case 2
                    spatialfilterstr = '4-NN variance filtering';
                case 3
                    spatialfilterstr = '27-NN variance filtering';
                case 4
                    spatialfilterstr = 'Gaussian smoothing';
                    spatialfilterstr = [spatialfilterstr ' [ ' nk_ConcatParamstr(params.SPATIAL.cubefwhm) ' ] ']; 
                case 5
                    spatialfilterstr = 'Resampling';
                    spatialfilterstr = [spatialfilterstr ' [ ' nk_ConcatParamstr(params.SPATIAL.cubevoxres) ' ] ']; 
            end
        end
        if isfield(params,'ACTPARAM') 
            
            for i=1:numel(params.ACTPARAM)
                
                if ~isempty(params.ACTPARAM{i})
                    
                    switch params.ACTPARAM{i}.cmd
                        
                        case 'labelimpute'
                            switch params.ACTPARAM{i}.LABELIMPUTE.method
                                case 'singlemean'
                                    impstr = 'single-subject median replacement';
                                case 'mean'
                                    impstr = 'Feature-wise mean replacement';
                                case 'dist'
                                    impstr = sprintf('median of %g nearest neighbors (%s)',params.ACTPARAM{i}.LABELIMPUTE.k, 'Manhattan');
                                case 'dist2'
                                    impstr = sprintf('median of %g nearest neighbors (%s)',params.ACTPARAM{i}.LABELIMPUTE.k, 'Euclidean');
                                case 'seuclidean'
                                    impstr = sprintf('median of %g nearest neighbors (%s)',params.ACTPARAM{i}.LABELIMPUTE.k, 'Seuclidean');
                                case 'cosine'
                                    impstr = sprintf('median of %g nearest neighbors (%s)',params.ACTPARAM{i}.LABELIMPUTE.k, 'Cosine similarity');
                                case 'mahalanobis'
                                    impstr = sprintf('median of %g nearest neighbors (%s)',params.ACTPARAM{i}.LABELIMPUTE.k, 'Mahalanobis');
                                case 'jaccard'
                                    impstr = sprintf('median of %g nearest neighbors (%s)',params.ACTPARAM{i}.LABELIMPUTE.k, 'Jaccard');
                                case 'hamming'
                                    impstr = sprintf('median of %g nearest neighbors (%s)',params.ACTPARAM{i}.IMPUTE.k, 'Hamming');
                            end
                            preprocact{i} = sprintf('Label Imputation: %s', impstr);
                            
                        case 'impute'
                            b = find(params.ACTPARAM{i}.IMPUTE.blockind);
                            switch params.ACTPARAM{i}.IMPUTE.method
                                case 'singlemean'
                                    impstr = 'single-subject median replacement';
                                case 'mean'
                                    impstr = 'Feature-wise mean replacement';
                                case 'dist'
                                    impstr = sprintf('median of %g nearest neighbors (%s)',params.ACTPARAM{i}.IMPUTE.k, 'Manhattan');
                                case 'dist2'
                                    impstr = sprintf('median of %g nearest neighbors (%s)',params.ACTPARAM{i}.IMPUTE.k, 'Euclidean');
                                case 'seuclidean'
                                    impstr = sprintf('median of %g nearest neighbors (%s)',params.ACTPARAM{i}.IMPUTE.k, 'Seuclidean');
                                case 'cosine'
                                    impstr = sprintf('median of %g nearest neighbors (%s)',params.ACTPARAM{i}.IMPUTE.k, 'Cosine similarity');
                                case 'mahalanobis'
                                    impstr = sprintf('median of %g nearest neighbors (%s)',params.ACTPARAM{i}.IMPUTE.k, 'Mahalanobis');
                                case 'jaccard'
                                    impstr = sprintf('median of %g nearest neighbors (%s)',params.ACTPARAM{i}.IMPUTE.k, 'Jaccard');
                                case 'hamming'
                                    impstr = sprintf('median of %g nearest neighbors (%s)',params.ACTPARAM{i}.IMPUTE.k, 'Hamming');
                                case 'hybrid'
                                    impstr = sprintf('median of %g nearest neighbors (%s)',params.ACTPARAM{i}.IMPUTE.k, 'Hybrid');    
                                otherwise
                                    impstr = sprintf('median of %g nearest neighbors (%s)',params.ACTPARAM{i}.IMPUTE.k, params.ACTPARAM{i}.IMPUTE.method); 
                            end
                            if isempty(b)
                                preprocact{i} = sprintf('Imputation in Matrix Block [ All features ]: %s', impstr);
                            else
                                preprocact{i} = sprintf('Imputation in Matrix Block [ %g columns, first: %g, last: %g ]: %s', numel(b), min(b), max(b), impstr);
                            end
                        case 'correctnuis'

                            if isfield(params.ACTPARAM{i},'METHOD') && params.ACTPARAM{i}.METHOD == 2
                                combatfl = 1;
                            else
                                combatfl = 0;
                            end
                                
                            if ~combatfl
                            
                                if ~isempty(params.ACTPARAM{i}.COVAR) 
                                    preprocact{i} = 'Partial correlations [ Covars:';
                                    covstr = ' ';
                                    for j=1:numel(params.ACTPARAM{i}.COVAR)
                                        covstr = [covstr res.covnames{params.ACTPARAM{i}.COVAR(j)} ', '];
                                    end
                                    preprocact{i} = [preprocact{i} covstr(1:end-2) ];
                                end
                                if isfield(params.ACTPARAM{i},'INTERCEPT')
                                    switch params.ACTPARAM{i}.INTERCEPT
                                        case 2
                                            preprocact{i} = [preprocact{i} ', intercept incl.' ];
                                        case 1
                                            preprocact{i} = [preprocact{i} ', intercept excl.' ];
                                    end
                                end
                                if isfield(params.ACTPARAM{i},'COVDIR')
                                    switch params.ACTPARAM{i}.COVDIR
                                        case 1
                                            preprocact{i} = [preprocact{i} ', remove effects' ];
                                        case 2
                                            preprocact{i} = [preprocact{i} ', increase effects' ];
                                    end
                                end
                            else
                                if ~isempty(params.ACTPARAM{i}.COVAR) 
                                      preprocact{i} = ['ComBat batch effect correction [ batch vector: ' res.covnames{params.ACTPARAM{i}.COVAR} ];
                                end
                                if params.ACTPARAM{i}.MCOVARUSE
                                    preprocact{i} = [ preprocact{i} ', retainment variable(s): ' strjoin( res.covnames(params.ACTPARAM{i}.MCOVAR),', ') ];
                                    if params.ACTPARAM{i}.MCOVARLABEL
                                        preprocact{i} = [ preprocact{i} ' + label' ];
                                    end
                                end
                                
                            end
                            preprocact{i} = [preprocact{i} ' ]'];
                            
                        case 'scale'

                            preprocact{i} = 'Scale';
                            if isfield(params.ACTPARAM{i}.SCALE,'AcMatFl')
                                if params.ACTPARAM{i}.SCALE.AcMatFl
                                    preprocact{i} = sprintf('%s across matrix', preprocact{i});
                                else
                                    preprocact{i} = sprintf('%s featurewise', preprocact{i});
                                end
                            end
                            if isfield(params.ACTPARAM{i}.SCALE,'ZeroOne')
                                if params.ACTPARAM{i}.SCALE.ZeroOne == 1
                                    preprocact{i} = sprintf('%s [ from 0 to 1 ]', preprocact{i});
                                else
                                    preprocact{i} = sprintf('%s [ from -1 to 1 ]', preprocact{i});
                                end
                            end
                            if params.ACTPARAM{i}.SCALE.zerooutflag == 1
                                preprocact{i} = sprintf('%s, zero-out completely non-finite features', preprocact{i});
                            end
                            
                        case 'normalize'
                            
                            if ~isempty(params.ACTPARAM{i}.IND) && ~strcmp(params.ACTPARAM{i}.IND,'NA')
                                preprocact{i} = 'Normalize data to group means (Group indices:';
                                covstr = ' ';
                                for j=1:numel(params.ACTPARAM{i}.IND)
                                    covstr = [covstr res.covnames{params.ACTPARAM{i}.IND(j)} ', '];
                                end
                                preprocact{i} = [preprocact{i} covstr(1:end-2) ')'];
                            else
                                preprocact{i} = 'Normalize data to group means (UNDEFINED)';
                            end
                            if params.ACTPARAM{i}.zerooutflag == 1
                                preprocact{i} = sprintf('%s, zero-out completely non-finite features', preprocact{i});
                            end
                        
                        case 'remmeandiff'
                            
                            preprocact{i} = 'Correct group offset from global mean => ';
                            if isfield(params.ACTPARAM{i},'sIND') && ...
                                    ~isempty(params.ACTPARAM{i}.sIND) && ...
                                    ~strcmp(params.ACTPARAM{i}.sIND,'NA')
                                preprocact{i} = [ preprocact{i} 'Source groups: ' ];
                                covstrS = ' ';
                                for j=1:numel(params.ACTPARAM{i}.sIND)
                                    covstrS = [covstrS res.covnames{params.ACTPARAM{i}.sIND(j)} ', '];
                                end
                                preprocact{i} = [preprocact{i} covstrS(1:end-2)];
                            end
                            if isfield(params.ACTPARAM{i},'dIND') && ...
                                    ~isempty(params.ACTPARAM{i}.dIND) && ...
                                    ~strcmp(params.ACTPARAM{i}.dIND,'NA')
                                preprocact{i} = [ preprocact{i} ', Destination groups: ' ];
                                covstrD = ' ';
                                for j=1:numel(params.ACTPARAM{i}.dIND)
                                    covstrD = [covstrD res.covnames{params.ACTPARAM{i}.dIND(j)} ', '];
                                end
                                preprocact{i} = [preprocact{i} covstrD(1:end-2)];
                            end
                            
                        case 'discretize'

                            preprocact{i} = ['Discretize (' ...
                                'Start at: ' num2str(params.ACTPARAM{i}.DISCRET.binstart) ...
                                ', Stepwidth: ' num2str(params.ACTPARAM{i}.DISCRET.binsteps) ...
                                ', Stop at: ' num2str(params.ACTPARAM{i}.DISCRET.binstop) ')' ];

                        case 'symbolize'

                            preprocact{i} = ['Symbolize (' ...
                                'Min. # bins:' num2str(params.ACTPARAM{i}.SYMBOL.symMinBin) ...
                                ', Max. # count:' num2str(params.ACTPARAM{i}.SYMBOL.symMaxBin) ...
                                ', Sequence length: ' num2str(params.ACTPARAM{i}.SYMBOL.symSeqLength) ...
                                ', #SDs: ' num2str(params.ACTPARAM{i}.SYMBOL.symStdNum) ')' ];

                        case 'reducedim'

                            if strcmp(params.ACTPARAM{i}.DR.RedMode,'PLS')
                                if isfield(params.ACTPARAM{i}.DR.PLS,'algostr')
                                    PLSalgo = sprintf('PLS (%s)',params.ACTPARAM{i}.DR.PLS.algostr);
                                else
                                    PLSalgo = 'PLS';
                                end
                                if params.ACTPARAM{i}.DR.PLS.uselabel==1
                                    RedMode = [ PLSalgo ': NM label' ];
                                else
                                    RedMode = [ PLSalgo ': ' num2str(size(params.ACTPARAM{i}.DR.PLS.V,2)) ' behavioral variables' ];
                                end
                                
                            else
                                RedMode = params.ACTPARAM{i}.DR.RedMode;
                                if isfield(params.ACTPARAM{i},'PX') && ~isempty(params.ACTPARAM{i}.PX)
                                    nParams = numel(params.ACTPARAM{i}.PX.Px); 
                                    paramstr = sprintf('[ %s: %s', ...
                                        params.ACTPARAM{i}.PX.Px(1).Params_desc, ...
                                        nk_ConcatParamstr(params.ACTPARAM{i}.PX.Px(1).Params));
                                    for j = 2:nParams
                                        paramstr = sprintf('%s: %s', ...
                                            params.ACTPARAM{i}.PX.Px(j).Params_desc, nk_ConcatParamstr(params.ACTPARAM{i}.PX.Px(j).Params));
                                    end
                                    RedMode = sprintf('%s, %s ]', RedMode, paramstr);
                                end
                            end
                            preprocact{i} = ['Reduce dimensionality (' RedMode ')'];
                            
                        case 'extdim'
                            preprocact{i} = 'Extract subspaces';
                            if isfield(params.ACTPARAM{i},'PX') && ~isempty(params.ACTPARAM{i}.PX)
                                nParams = numel(params.ACTPARAM{i}.PX.Px); 
                                paramstr = sprintf('[ %s: %s', ...
                                    params.ACTPARAM{i}.PX.Px(1).Params_desc, ...
                                    nk_ConcatParamstr(params.ACTPARAM{i}.PX.Px(1).Params));
                                for j = 2:nParams
                                    paramstr = sprintf('%s, %s: %s', paramstr, ...
                                        params.ACTPARAM{i}.PX.Px(j).Params_desc, nk_ConcatParamstr(params.ACTPARAM{i}.PX.Px(j).Params));
                                end
                                preprocact{i} = sprintf('%s %s ]', preprocact{i}, paramstr);
                            end
                            
                        case 'standardize'
                            
                            if isfield(params.ACTPARAM{i},'METHOD')
                                METHOD = params.ACTPARAM{i}.METHOD;
                            else
                                METHOD = 'Standardization';
                            end
                            preprocact{i} = METHOD;
                            if ~isempty(params.ACTPARAM{i}.WINSOPT)
                                preprocact{i} = sprintf('%s & winsorization at +/- %s', ...
                                    preprocact{i}, nk_ConcatParamstr(params.ACTPARAM{i}.WINSOPT));
                            end
                            if isfield(params.ACTPARAM{i},'sIND') && ~isempty(params.ACTPARAM{i}.sIND)
                                preprocact{i} = sprintf('%s, computed in subgroup (%g cases)', ...
                                    preprocact{i}, sum(params.ACTPARAM{i}.sIND));
                            end
                            if isfield(params.ACTPARAM{i},'dIND') && ~isempty(params.ACTPARAM{i}.dIND)
                                preprocact{i} = sprintf('%s, applied to subgroup (%g cases)', ...
                                    preprocact{i}, sum(params.ACTPARAM{i}.dIND));
                            end
                            if params.ACTPARAM{i}.zerooutflag == 1
                                preprocact{i} = sprintf('%s, zero-out completely non-finite features', preprocact{i});
                            end
                            
                        case 'unitnormalize'
                            
                            preprocact{i} = 'Normalize to unit vector';
                            switch params.ACTPARAM{i}.METHOD
                                case 1
                                    preprocact{i} = [ preprocact{i} ' (L1-norm)' ];
                                case 2
                                    preprocact{i} = [ preprocact{i} ' (L2-norm)' ];
                            end
                             if params.ACTPARAM{i}.zerooutflag == 1
                                preprocact{i} = sprintf('%s, zero-out completely non-finite features', preprocact{i});
                            end

                        case 'elimzero'
                            
                            preprocact{i} = 'Prune non-informative columns from matrix [';
                            if params.ACTPARAM{i}.PRUNE.zero == 1
                                preprocact{i} = sprintf('%s Zero Var', preprocact{i});
                            end
                            if params.ACTPARAM{i}.PRUNE.nan == 1
                                preprocact{i} = sprintf('%s, Nan', preprocact{i});
                            end    
                            if params.ACTPARAM{i}.PRUNE.inf == 1
                                preprocact{i} = sprintf('%s, Inf', preprocact{i});
                            end
                            if isfield(params.ACTPARAM{i}.PRUNE,'perc') && ~isempty(params.ACTPARAM{i}.PRUNE.perc)
                                preprocact{i} = sprintf('%s, Single-Value [ %s ]', preprocact{i}, nk_ConcatParamstr(params.ACTPARAM{i}.PRUNE.perc));
                            end
                            preprocact{i} = sprintf('%s ]', preprocact{i});
                             
                        case 'rankfeat'
                            
                            switch params.ACTPARAM{i}.RANK.algostr 
                                case 'pearson'
                                    if params.ACTPARAM{i}.RANK.Pearson == 1, algostr = 'Pearson'; else, algostr = 'Spearman'; end
                                case 'extern'
                                    algostr = 'External ranking';
                                case 'extern_fscore'
                                    algostr = 'External combined with f-score ranking';
                                otherwise
                                    algostr = params.ACTPARAM{i}.RANK.algostr;
                            end
                            
                            switch params.ACTPARAM{i}.RANK.weightmethod
                                case 1
                                    preprocact{i} = sprintf('Rank up features using %s ', algostr );
                                case 2
                                    preprocact{i} = sprintf('Rank down features using %s ', algostr );
                            end
                            
                            if ~any(strcmp(params.ACTPARAM{i}.RANK.algostr,{'extern','pls'}))
                                switch params.ACTPARAM{i}.RANK.ranktype
                                    case 1
                                         preprocact{i} = sprintf( '%s => target label', preprocact{i} );
                                    case 2
                                         preprocact{i} = sprintf( '%s => categorical label: %s', preprocact{i}, params.ACTPARAM{i}.RANK.labeldesc );
                                    case 3
                                         preprocact{i} = sprintf( '%s => continuous label: %s', preprocact{i}, params.ACTPARAM{i}.RANK.labeldesc);
                                end
                            end
                            
                        case 'extfeat'
                            if isfield(params.ACTPARAM{i}.W_ACT,'softflag')
                                switch params.ACTPARAM{i}.W_ACT.softflag
                                    case 1
                                        threshstr = sprintf('exp. multiplier(s): %s', nk_ConcatParamstr(params.ACTPARAM{i}.W_ACT.exponent));
                                        threshstr = [ 'Soft feature selection (' threshstr ')' ];
                                    case 2
                                        threshstr = sprintf('threshold(s): %s', nk_ConcatParamstr(params.ACTPARAM{i}.W_ACT.threshvec));
                                        if params.ACTPARAM{i}.W_ACT.clustflag == 1, 
                                            threshstr = [ 'clusterized ' threshstr ]; 
                                        end
                                        threshstr = [ 'Hard feature selection (' threshstr ')' ];

                                end
                            else
                                threshstr = 'undefined';
                            end
                            preprocact{i} = threshstr;
                            
                        case 'remvarcomp'
                            if isfield(params.ACTPARAM{i}.REMVARCOMP,'G')
                                [p,q] = size(params.ACTPARAM{i}.REMVARCOMP.G);
                                if q>1
                                    REMVARCOMP_G_str = sprintf('Target: %g-by-%g Matrix',p,q);
                                else
                                    REMVARCOMP_G_str = 'Target: Vector';
                                end
                            else
                                REMVARCOMP_G_str = 'No covariates specified';
                            end
                            switch params.ACTPARAM{i}.REMVARCOMP.corrmeth
                                case 1
                                    REMVARCOMP_corrmeth_str = 'Correlation Method: Pearson';
                                case 2
                                    REMVARCOMP_corrmeth_str = 'Correlation Method: Spearman';
                                case 3
                                    REMVARCOMP_corrmeth_str = 'Correlation Method: ANOVA';
                            end
                            if isfield(params.ACTPARAM{i}.REMVARCOMP,'recon')
                                switch params.ACTPARAM{i}.REMVARCOMP.recon
                                    case 1
                                        REMVARCOMP_recon_str = 'yes';
                                    case 2
                                        REMVARCOMP_recon_str = 'no';
                                end
                            else
                                REMVARCOMP_recon_str = 'undefined';
                            end
                            if isfield(params.ACTPARAM{i}.REMVARCOMP,'varop')
                                REMVARCOMP_varop_str = params.ACTPARAM{i}.REMVARCOMP.varop;
                            else
                                REMVARCOMP_varop_str = 'undefined';
                            end
                            REMVARCOMP_corrthresh_str = nk_ConcatParamstr(params.ACTPARAM{i}.REMVARCOMP.corrthresh);
                            preprocact{i} = sprintf('Variance extraction [ %s, %s, Correlation cutoff(s): %s, Operator: %s, Back-projection: %s ]', ...
                                REMVARCOMP_G_str, REMVARCOMP_corrmeth_str , REMVARCOMP_corrthresh_str, REMVARCOMP_varop_str, REMVARCOMP_recon_str);
                        case 'devmap'
                            
                            if ~isempty(params.ACTPARAM{i}.DEVMAP.glabel)
                                grpstr = sprintf(', computation in %g cases', sum(params.ACTPARAM{i}.DEVMAP.glabel));
                            else
                                grpstr = '';
                            end
                            if size(params.ACTPARAM{i}.DEVMAP.covmat,2)>1
                                covsizeextr = sprintf('%g covariates, %g components', size(params.ACTPARAM{i}.DEVMAP.covmat,2), params.ACTPARAM{i}.DEVMAP.ncomp);   
                            else
                                covsizeextr = sprintf('1 covariate');
                            end
                            preprocact{i} = sprintf('Deviation-based weighting [ %s: %s%s ]', params.ACTPARAM{i}.DEVMAP.algostr, covsizeextr, grpstr);
                                
                    end
                else
                    preprocact{i} = 'undefined';
                end
            end
        end
        vargout.PREPROC.groupmode = groupmode;
        vargout.PREPROC.targetscaling = targetscalestr;
        vargout.PREPROC.labelimpute = labelimputestr;
        vargout.PREPROC.spatialfiltering = spatialfilterstr;
        vargout.PREPROC.preprocact = preprocact;
        
    case 'FeatFlt'
        vargout.FeatFltFlag = 'no';
        if isfield(params.Filter,'EnsembleStrategy')
         switch params.Filter.EnsembleStrategy.Metric
            case 1
                vargout.FeatMetric = 'Predicted Targets';
            case 2
                vargout.FeatMetric = 'ML Algorithm scores';
         end 
        else
            vargout.FeatMetric = 'NA';
        end
        if params.Filter.flag
            vargout.FeatFltFlag = 'yes';
            if params.Filter.binmode
                vargout.FeatFltBinmode = 'yes';
            else
                vargout.FeatFltBinmode = 'no';
            end
            if params.Filter.SubSpaceFlag
                switch params.Filter.SubSpaceStrategy
                    case 1 
                        SubSpaceStrat = 'Winner takes all';
                    case 2
                        SubSpaceStrat = sprintf('%g%% range from optimum',params.Filter.SubSpaceCrit);
                    case 3
                        SubSpaceStrat = sprintf('%g%%-percentile',params.Filter.SubSpaceCrit);
                    case 4
                        SubSpaceStrat = 'All subspaces';
                end
                vargout.FeatFltSubSpaces = sprintf('Feature Subspace Optimization enabled ( %s )', SubSpaceStrat);
            else
                vargout.FeatFltSubSpaces = 'Ranking Optimization enabled';
                vargout.FeatFltThresh = sprintf('Threshold: %1.2f',params.Filter.RankThresh);
            end
            switch params.Filter.SubSpaceStepping
                case 0
                    vargout.FeatFltStepping = 'Every feature will be evaluated';
                otherwise
                    vargout.FeatFltStepping = sprintf('Blocks with %1.0f%% of features', params.Filter.SubSpaceStepping);
            end
            switch params.Filter.CostFun
                case 1
                    vargout.FeatPop  = 'CV1 Training Data';
                case 2
                    vargout.FeatPop  = 'CV1 Test Data';
                case 3
                    vargout.FeatPop  = 'CV1 Training & Test Data';
            end
            switch params.Filter.type
                case 0
                    vargout.FeatFlt = 'none';
                case 1
                    vargout.FeatFlt = Get_FEASTDescription(params.Filter);
                    vargout.FeatFlt = ['FEAST (' vargout.FeatFlt ')'];
                case 2
                    vargout.FeatFlt = 'MRMR';
                case 3
                    if isfield(params.Filter,'Pearson'),
                        switch params.Filter.Pearson
                            case 1
                                vargout.FeatFlt = 'Pearson';
                            case 2
                                vargout.FeatFlt = 'Spearman';
                        end
                    else
                        vargout.FeatFlt = 'Pearson';
                    end
                case 4
                    vargout.FeatFlt = 'Simba';
                case 5
                    vargout.FeatFlt = Get_GflipDescription(params.Filter);
                    vargout.FeatFlt = ['Gflip [ ' vargout.FeatFlt ']'];
                case 6
                    vargout.FeatFlt = 'AMS';
                case 7
                    vargout.FeatFlt = Get_IMReliefDescription(params.Filter);
                    vargout.FeatFlt = ['IMRelief [ ' vargout.FeatFlt ']'];
                case 9
                    vargout.FeatFlt = 'Incr. card.';
                case 10
                    vargout.FeatFlt = Get_RGSDescription(params.Filter);
                    vargout.FeatFlt = ['RGS [ ' vargout.FeatFlt ']'];
                case 11
                    vargout.FeatFlt = Get_RSSDescription(params.Filter);
                case 12
                    vargout.FeatFlt = ['Relief [ k = ' num2str(params.Filter.Relief.k) ' ]' ];
                case 13
                    vargout.FeatFlt = 'FScore';
                case 14
                    vargout.FeatFlt = 'Bhattacharya distance';
            end
                
            if params.Filter.SubSpaceStrategy > 1 && params.Filter.SubSpaceFlag
                
                switch params.Filter.EnsembleStrategy.type
                    case 0
                        vargout.FeatFltEnsStrat = 'Aggregated Ensemble';

                    case 9
                        vargout.FeatFltEnsStrat = ['Probabilistic Feature Extraction (agreement: ' ...
                            num2str(params.Filter.EnsembleStrategy.Perc) ...
                            ', min # of features: ' ...
                            num2str(params.Filter.EnsembleStrategy.MinNum) ')'];

                    otherwise 

                        switch params.Filter.EnsembleStrategy.type
                            case 1
                                vargout.FeatFltEnsStrat = 'Recursive Ensemble Thinning (entropy: ';
                            case 2
                                vargout.FeatFltEnsStrat = 'Recursive Ensemble Thinning (performance: ';
                            case 3
                                vargout.FeatFltEnsStrat = 'Recursive Ensemble Thinning (diversity: ';
                            case 4
                                vargout.FeatFltEnsStrat = 'Recursive Ensemble Thinning (bias & variance: ';
                            case 5
                                vargout.FeatFltEnsStrat = 'Forward Ensemble Construction (entropy: ';
                            case 6
                                vargout.FeatFltEnsStrat = 'Forward Ensemble Construction (performance: ';
                            case 7
                                vargout.FeatFltEnsStrat = 'Forward Ensemble Construction (diversity: ';
                            case 8
                                vargout.FeatFltEnsStrat = 'Forward Ensemble Construction (bias & variance: ';
                        end

                        vargout.FeatFltEnsStrat = [vargout.FeatFltEnsStrat vargout.FeatPop];

                end 
                vargout.FeatFltEnsStrat = [vargout.FeatFltEnsStrat ', Metric: ' vargout.FeatMetric ')'];
                
            else
                vargout.FeatFltEnsStrat = 'No subspace ensemble construction (winner takes it all)';
                
            end
            %if isfield(RFE.Filter,'SubSpaceCrit')
            vargout.FilterMode = [vargout.FeatFlt ' => ' vargout.FeatFltSubSpaces];
            vargout.FilterMethod = ['Method: ' vargout.FeatFltEnsStrat];
        else
            vargout.FilterMode = 'none';
            vargout.FilterMethod = '';
        end
               
    case 'EnsType'
        
        switch params.CostFun
            case 1
                vargout.CostType = 'argmax(CV=>training data)';
            case 2
                vargout.CostType = 'argmax(CV=>test data)';
            case 3
                vargout.CostType = 'argmax(CV=>training & test data)';
        end
        switch params.SubSpaceStrategy
            case 1
                vargout.SubSpaceStrat = 'Winner subspace takes it all';
            case {2,3}
                 switch params.SubSpaceStrategy
                    case 2
                        vargout.SubSpaceStrat = sprintf('Subspace ensemble in a range of %g from the optimum', params.SubSpaceCrit);
                    case 3
                        vargout.SubSpaceStrat = sprintf('Subspace ensemble at the %g-percentile ', params.SubSpaceCrit);
                 end
            case 4
                vargout.SubSpaceStrat = 'Ensemble of all subspaces';
        end
        if isfield(params.EnsembleStrategy,'Perc')
            vargout.EnsPerc = sprintf('%g',params.EnsembleStrategy.Perc);
        else
            vargout.EnsPerc = 'NA';
        end
        if isfield(params.EnsembleStrategy,'MinNum')
            vargout.EnsMinNum = sprintf('%g',params.EnsembleStrategy.MinNum);
        else
            vargout.EnsMinNum = 'NA';
        end
        switch params.EnsembleStrategy.DataType
            case 1
                vargout.FeatPop     = 'CV1 Training Data';
            case 2
                vargout.FeatPop     = 'CV1 Test Data';
            case 3
                vargout.FeatPop     = 'CV1 Training & Test Data';
        end
        switch params.EnsembleStrategy.Metric
            case 1
                vargout.FeatMetric = 'Predicted Targets';
            case 2
                vargout.FeatMetric = 'ML Algorithm scores';
        end
        switch params.EnsembleStrategy.Weighting
            case 1
                vargout.EnsWeighting = 'Enabled'; 
            case 0
                vargout.EnsWeighting = 'Disabled';
        end
        switch params.EnsembleStrategy.ConstructMode
            case 0
                vargout.EnsConMode = 'Aggregated (Bagged) Ensemble';
            case 1
                vargout.EnsConMode = 'Greedy Backward Base Learner Elimination';
            case 2
                vargout.EnsConMode = 'Greedy Forward Ensemble Construction';
            case 4
                vargout.EnsConMode = 'Ensemble-based Probabilistic Feature Extraction';
            case 5
                vargout.EnsConMode = 'Adaboost (NOT TESTED)';
        end
        switch params.EnsembleStrategy.DivCrit
            case 0
                vargout.EnsDivCrit = 'NA';
            case 1
                vargout.EnsDivCrit = 'Entropy'; 
            case 2
                vargout.EnsDivCrit = 'Entropy & Performance'; 
            case 3
                vargout.EnsDivCrit = 'Kappa Diversity'; 
            case 4
                vargout.EnsDivCrit = 'Bias-Variance decomposition'; 
            case 5
        end
        switch params.EnsembleStrategy.type
            case 0
                vargout.EnsStrat = 'Aggregated Ensemble';
            case 9
                vargout.EnsStrat = ['Probabilistic Feature Extraction (agreement: ' ...
                    vargout.EnsPerc ', min # of features: ' vargout.EnsMinNum ')'];
            otherwise 
                vargout.EnsStrat = [vargout.EnsConMode ' (' vargout.EnsDivCrit ')' ];
                vargout.EnsStrat = [vargout.EnsStrat ', ' vargout.FeatPop];
        end 
        vargout.EnsStrat = [vargout.EnsStrat ', Metric: ' vargout.FeatMetric ')'];
        
    case 'FeatWrap'
        
        if isfield(params.Wrapper,'optflag')
            if params.Wrapper.optflag == 1
                vargout.WrapperOptFlag = 'Wrapper at parameter optimum only';
            else
                vargout.WrapperOptFlag = 'Wrapper at all parameter combinations';
            end
        else
            vargout.WrapperOptFlag = 'undefined';
        end
        
        if params.Wrapper.flag || params.Wrapper.optflag == 1
                    
            if params.Wrapper.PFE.flag
                vargout.WrapperPFE = 'Cross-CV1 PFE on';
                if ~isfield(params.Wrapper.PFE,'Mode'), params.Wrapper.PFE.Mode = 3; end
                switch params.Wrapper.PFE.Mode
                    case 1
                        vargout.WrapperPFE = sprintf('%s [ thresh: %g% (-%g), Min. No. of feats: %s]',vargout.WrapperPFE, ...
                                params.Wrapper.PFE.Perc, ...
                                params.Wrapper.PFE.TolWin, ...
                                params.Wrapper.PFE.MinNum);
                    case 2
                        vargout.WrapperPFE = sprintf('%s [ %g% top feats ]',vargout.WrapperPFE,params.Wrapper.PFE.Perc);
                    case 3
                        vargout.WrapperPFE = sprintf('%s [ %g%% top feats ]',vargout.WrapperPFE,params.Wrapper.PFE.Perc);
                end
            else
                vargout.WrapperPFE = sprintf('Cross-CV1 PFE off');
            end
            
            switch params.Wrapper.datamode
                case 1
                    vargout.WrapperDataMode = 'CV1 training data';
                case 2
                    vargout.WrapperDataMode = 'CV1 test data';
                case 3
                    vargout.WrapperDataMode = 'CV1 training & test data';
            end
            if isfield(params.Wrapper,'GreedySearch')
                switch params.Wrapper.GreedySearch.FeatStepPerc
                    case 0
                        vargout.WrapperFeatStepPerc = 'Every feature';
                    otherwise
                        vargout.WrapperFeatStepPerc = sprintf('%g%% of features',params.Wrapper.GreedySearch.FeatStepPerc);
                end

                switch params.Wrapper.GreedySearch.EarlyStop.Thresh
                    case 0
                        vargout.WrapperEarlyStop = sprintf('EarlyStop disabled');                
                    otherwise 
                        switch params.Wrapper.GreedySearch.EarlyStop.Perc
                            case 1
                                vargout.WrapperEarlyStop = sprintf('Stop at k=%g%% of features',params.Wrapper.GreedySearch.EarlyStop.Thresh);
                            case 2
                                vargout.WrapperEarlyStop = sprintf('Stop at k=%g features',params.Wrapper.GreedySearch.EarlyStop.Thresh);
                        end
                end
            end
            
            switch params.Wrapper.type
                case 1
                    vargout.WrapperMode = 'Greedy feature selection';
                    vargout.WrapperStr = sprintf('%s (Use %s; %s; %s)',vargout.WrapperMode, vargout.WrapperDataMode, vargout.WrapperEarlyStop, vargout.WrapperPFE);
                case 2
                    vargout.WrapperMode = 'Simulated annealing';
                    vargout.WrapperStr = sprintf('%s',vargout.WrapperMode);
            end
            
            if params.Wrapper.SubSpaceStrategy > 1
                switch params.Wrapper.EnsembleStrategy.type
                    case 0
                        vargout.FeatWrapEnsStrat = 'No subspace evaluation';
                    case {1, 2, 3, 5, 6, 7, 9}

                        switch params.Wrapper.EnsembleStrategy.type
                            case 1
                                vargout.FeatWrapEnsStrat = 'Recursive Ensemble Thinning (entropy: ';
                            case 2
                                vargout.FeatWrapEnsStrat = 'Recursive Ensemble Thinning (performance: ';
                            case 3
                                vargout.FeatWrapEnsStrat = 'Recursive Ensemble Thinning (diversity: ';
                            case 5
                                vargout.FeatWrapEnsStrat = 'Forward Ensemble Construction (entropy: ';
                            case 6
                                vargout.FeatWrapEnsStrat = 'Forward Ensemble Construction (performance: ';
                            case 7
                                vargout.FeatWrapEnsStrat = 'Forward Ensemble Construction (diversity: ';
                            case 9 
                                vargout.FeatWrapEnsStrat = 'No Ensemble Method used (Optimize using ';
                        end

                        switch params.Wrapper.EnsembleStrategy.DataType
                            case 1
                                vargout.FeatWrapEnsStrat = [vargout.FeatWrapEnsStrat 'Training Data)'];
                            case 2
                                vargout.FeatWrapEnsStrat = [vargout.FeatWrapEnsStrat 'CV1 Test Data)'];
                        end

                    case 4
                        vargout.FeatWrapEnsStrat = ['Probabilistic Feature Extraction (agreement: ' ...
                            num2str(params.Wrapper.EnsembleStrategy.Perc) ...
                            ', min # of features: ' ...
                            num2str(params.Wrapper.EnsembleStrategy.MinNum) ')'];

                end 

                switch params.Wrapper.EnsembleStrategy.Metric
                        case 1
                            vargout.FeatWrapEnsStrat = [vargout.FeatWrapEnsStrat ', Metric: Predicted targets)'];
                        case 2
                            vargout.FeatWrapEnsStrat = [vargout.FeatWrapEnsStrat ', Metric: Dec. Values / Prob.)'];
                end
                
            else
                vargout.FeatWrapEnsStrat = 'No subspace ensemble construction (winner takes it all)';
            end
             vargout.WrapperMethod = ['Subspace method: ' vargout.FeatWrapEnsStrat];
        else
            vargout.WrapperMode = 'none';
            vargout.WrapperStr = 'none';
            vargout.WrapperMethod = '';
        end
        
	case 'GridParam'
        if isfield(params.SVM,'GridParam')
            switch params.SVM.GridParam
                case 1
                    vargout.GridParam = 'Accuracy';
                case 2
                    vargout.GridParam = 'TPR';
                case 3
                    vargout.GridParam = 'Sensitivity';
                case 4
                    vargout.GridParam = 'FPR';
                case 5
                    vargout.GridParam = 'PPV';
                case 6
                    vargout.GridParam = 'MCC';
                case 7
                    vargout.GridParam = 'AUC';
                case 9
                    vargout.GridParam = 'MSE';
                case 10
                    vargout.GridParam = 'SCC (R^2)';
                case 11
                    vargout.GridParam = 'NRMSD';
                case 12
                    vargout.GridParam = 'RMSD';
                case 13
                    vargout.GridParam = 'Gmean';
                case 14
                    vargout.GridParam = 'BAC';
                case 15
                    vargout.GridParam = 'F-score';
                case 16
                    vargout.GridParam = 'CC';
                case 17
                    vargout.GridParam = 'Sens*Spec';
                case 18
                    vargout.GridParam = 'MAE'; 
                case 19
                    vargout.GridParam = 'PSI';
                case 20
                    vargout.GridParam = 'NNP';
            end
        else
            vargout.GridParam = 'NA';
        end
        
    case 'ParamComb'
        
        PX_preML = []; preML = []; ML = []; preML_nCombs = 0; ML_nCombs = 0;
        if isfield(params,'PREPROC')
            PX_preML = nk_ReturnParamChain(params.PREPROC(varind), true); nP = numel(PX_preML);
        end
        PX_ML = nk_ReturnParamChain(params.GRD, true);
        cnt = 1;
        for n=1:nP
            if nP>1, PXPREML = PX_preML(n); else PXPREML = PX_preML; end
            if isfield(PXPREML,'Params') && ~isempty(PXPREML.Params)
                for i=1:numel(PXPREML.Params)
                    if nP>1, ModalStr = sprintf(' [ Modality #%g ]', n); else ModalStr = ''; end
                    preML{cnt} = sprintf('[ %s%s => %s ]', PXPREML.Params_desc{i}, ModalStr, nk_ConcatParamstr(PXPREML.Params{i}));
                    cnt = cnt+1;
                end
                preML_nCombsn = size(PX_preML(n).opt,1);
            else
                preML{cnt} = sprintf('none');
                cnt = cnt+1;
                preML_nCombsn = 0; 
            end
            if preML_nCombs == 0, 
                preML_nCombs = preML_nCombsn;
            else
                preML_nCombs = preML_nCombs * preML_nCombsn;
            end
        end
        vargout.preML_nCombs = preML_nCombs;
        if isfield(PX_ML,'Params') && ~isempty(PX_ML.Params)
            for i=1:numel(PX_ML.Params)
                ML{i} = sprintf('[ %s => %s ]', PX_ML.Params_desc{i},nk_ConcatParamstr(PX_ML.Params{i}));
            end
            ML_nCombsn = size(PX_ML.opt,1); 
        else
            ML{1} = sprintf('none');
            ML_nCombsn = 0; 
        end  
        if ML_nCombs == 0, 
            ML_nCombs = ML_nCombsn;
        else
            ML_nCombs = ML_nCombs * ML_nCombsn;
        end
        vargout.ML_nCombs = ML_nCombs;
        vargout.preML = preML; vargout.ML = ML;
        
    case 'SVMprog'
        switch params.SVM.prog
            case 'LIBSVM'
                switch params.SVM.LIBSVM.LIBSVMver
                    case 0
                        vargout.prog = 'LIBSVM 3.1.2 with instance weighting support';
                    case 1
                        vargout.prog = 'LIBSVM 2.9.1';
                    case 2
                        vargout.prog = 'LIBSVM 2.89';
                    case 3
                        vargout.prog = 'LIBSVM 2.89 PLUS';
                end
            case 'LIBLIN'
                vargout.prog = 'LIBLINEAR';    
            case 'SVMLIT'
                vargout.prog = 'SVMlight';
            case 'SVMPRF'
                vargout.prog = 'SVMperf';
            case 'MikRVM'
                vargout.prog = 'RVM (Mike Tipping)';
            case 'MKLRVM'
                vargout.prog = 'RVM (Psorakis & Damoulas)';
            case 'IMRELF'
                vargout.prog = 'IMRelief';
            case 'GLMFIT'
                vargout.prog = 'GLM';
            case 'LIKNON'
                vargout.prog = 'LIKNON dichotomizer';
            case 'kNNMEX'
                vargout.prog = 'kNN';
            case 'BLOREG'
                vargout.prog = 'Sparse Bayesian Logistic Regression';
            case 'LSTSVM'
                vargout.prog = 'Least-Squares Support Vector Machine';
            case 'MVTRVR'
                vargout.prog = 'Multinomial Relevance Vector Regression';
            case 'MVTRVM'
                vargout.prog = 'Multinomial Relevance Vector Classification';
            case 'MEXELM'
                vargout.prog = 'Extreme Learning Machine';
            case 'MSTOOL'
                vargout.prog = 'Marc Schmidt ML algorithm';
            case 'FAMRVR'
                vargout.prog = 'Fast Multivariate Relevance Vector Regression';
            case 'DECTRE'
                vargout.prog = 'Decision Tree';
            case 'RNDFOR'
                vargout.prog = 'Random Forest';
            case 'matLRN'
                vargout.prog = 'matLRN Library';
            case 'SEQOPT'
                vargout.prog = 'Sequence Optimizer';
            otherwise
                vargout.prog = params.SVM.prog;
        end
        
    case 'classifier'
        switch params.SVM.prog
            case 'IMRelief'
                vargout.classifier = params.SVM.imrelief.distance;
            case {'MikRVM','MKLRVM'}
                vargout.classifier = 'RVM';
            case {'LIBSVM'}
                classifier = params.SVM.LIBSVM.classifier;
                switch classifier
                    case 0
                        vargout.classifier = 'L1-Loss SVC';
                    case 1
                        switch params.SVM.LIBSVM.LIBSVMver
                            case {0,1,2}
                                vargout.classifier = 'nu-SVC';
                            case 3
                                vargout.classifier = 'L2-Loss SVC';
                        end
                    case 2
                        switch params.SVM.LIBSVM.LIBSVMver
                            case {0,1,2}
                                vargout.classifier = 'one-class SVM';
                            case 3
                                vargout.classifier = 'nu-SVC';
                        end
                    case 3
                        switch params.SVM.LIBSVM.LIBSVMver
                            case {0,1,2}
                                classifier = 'eps-SVR';
                                epsparam    = params.SVM.LIBSVM.Optimization.p;
                                vargout.classifier = [classifier '(' num2str(epsparam) ')'];
                            case 3
                                vargout.classifier = 'one-class SVM';
                        end

                    case 4
                        switch params.SVM.LIBSVM.LIBSVMver
                            case {0,1,2}
                                classifier = 'nu-SVR';
                                nuparam     = params.SVM.LIBSVM.Optimization.nu;
                                vargout.classifier = [classifier '(' num2str(nuparam) ')'];
                            case 3
                                classifier = 'eps-SVR';
                                epsparam    = params.SVM.LIBSVM.Optimization.p;
                                vargout.classifier = [classifier '(' num2str(epsparam) ')'];
                        end
                    case 5
                        vargout.classifier = 'nu-SVR';
                                nuparam     = params.SVM.LIBSVM.Optimization.nu;
                                vargout.classifier = [classifier '(' num2str(nuparam) ')'];
                    case 6
                        vargout.classifier = 'Support-Vector Domain Descriptor (L1-Loss SVC)';
                    case 7
                        vargout.classifier = 'Support-Vector Domain Descriptor (L2-Loss SVC)';
                    otherwise
                        vargout.classifier = ['Old NeuroMiner version: ' classifier];
                end
            case 'LIBLIN'
                if isfield(params.SVM,'LIBLIN')
                    switch params.SVM.LIBLIN.classifier
                        case 0
                            vargout.classifier = 'L2-regularized logistic regression (primal)';
                        case 1
                            vargout.classifier =  'L2-regularized L2-loss SVC (dual)';
                        case 2
                            vargout.classifier =  'L2-regularized L2-loss SVC (primal)';
                        case 3
                            vargout.classifier =  'L2-regularized L1-loss SVC (dual)';
                        case 5
                            vargout.classifier =  'L1-regularized L2-loss SVC';
                        case 6
                            vargout.classifier =  'L1-regularized logistic regression)';
                        case 7
                            vargout.classifier = 'L2-regularized logistic regression (dual)';
                        case 11
                            vargout.classifier = 'L2-regularized L2-loss SVR (primal)';
                        case 12
                            vargout.classifier = 'L2-regularized L2-loss SVR (dual)';
                        case 13
                            vargout.classifier = 'L2-regularized L1-loss SVR (dual)';
                    end
                    vargout.classifier = [ vargout.classifier ', Tolerance: ' num2str(params.SVM.LIBLIN.tolerance) ];                
                else
                    vargout.classifier = 'undefined';
                end
            case 'GLMFIT'
                vargout.classifier = 'Linear Regression Model';
            case 'LIKNON'
                vargout.classifier = 'Sparse Linear Hyperplane';
            case 'kNNMEX'
                vargout.classifier = 'kNN decision rule';
            case 'IMRELF'
                vargout.classifier = 'Iterative Local Linear Optimization';
            case 'BLOREG'
                vargout.classifier = 'blogreg';
            case 'LSTSVM'
                vargout.classifier = params.SVM.LSTSVM.type;
            case {'MVTRVR','MVTRVM'}
                if isfield(params.SVM,'MVTRVR')
                    vargout.classifier = ['RVM with ' num2str(params.SVM.MVTRVR.iter) ' iterations'];
                else
                    vargout.classifier = 'RVM (undefined parameters)';
                end
            case 'FAMRVR'
                vargout.classifier = ['Fast RVR with ' num2str(params.SVM.FAMRVR.iter) ' iters; tolerance = '  num2str(params.SVM.FAMRVR.tolerance)];
            case 'MEXELM'
                vargout.classifier = 'Extreme Learning Machine (ELM)';
            case 'DECTRE'
                vargout.classifier = 'MATLAB''s fitctree / fitrtree algorithm';
            case 'RNDFOR'
                vargout.classifier = 'MATLAB''s Random Forest algorithm';
            case 'matLRN'
                if isfield(params.SVM,'matLRN')
                    vargout.classifier = sprintf('ml_%s_%s', params.SVM.matLRN.learner.framework, char(params.SVM.matLRN.algo));
                else
                    vargout.classifier = 'matLearn setup not available';
                end
            case 'GLMNET'
                vargout.classifier = 'LASSO/Elastic-net regularized GLM';
            case 'GRDBST'
                vargout.classifier = 'Gradient Boosting';
            case 'ROBSVM'
                vargout.classifier = 'Robust LIBSVM';
            case 'SEQOPT'
                vargout.classifier = 'Sequence optimization based on uncertainty case propagation';
            case 'WBLCOX'
                vargout.classifier = 'Willbur-Cox proportional hazard regression';
        end
        
    case 'kernel'            
        vargout.kernel = params.SVM.kernel.kerndesc;
        
    case 'cv'
        vargout.cv = ['CV2: ' num2str(size(params.cv.TrainInd,1)) 'x' num2str(size(params.cv.TrainInd,2)) ...
                        ', CV1: ' num2str(size(params.cv.cvin{1,1}.TrainInd,1)) 'x' num2str(size(params.cv.cvin{1,1}.TrainInd,2))];
        vargout.cv2samples = size(params.cv.TrainInd,1) * size(params.cv.TrainInd,2) ;
        vargout.cv1samples = vargout.cv2samples * size(params.cv.cvin{1,1}.TrainInd,1) * size(params.cv.cvin{1,1}.TrainInd,2);
        
    case 'multiclass'
        vargout.multiclass = 'NA';
        
        if isfield(params,'MULTI') 
            if params.MULTI.flag
                vargout.multiclass = 'Enabled';

                if isfield(params.MULTI,'train') && params.MULTI.train
                    vargout.multiclass = [vargout.multiclass ', Multi-group optimization'];
                else
                    vargout.multiclass = [vargout.multiclass ', Binary optimization'];
                end
                if isfield(params.MULTI,'method')
                    switch params.MULTI.method
                        case 1 % One-Vs-One Max Wins
                            vargout.multiclass = [vargout.multiclass ', One-vs-One-Max-Wins'];
                            switch params.MULTI.decisiontype
                                case 1
                                    decisiontypestr = 'Sum';
                                case 2
                                    decisiontypestr = 'Mean';
                                case 3
                                    decisiontypestr = 'Product';
                                case 4
                                    decisiontypestr = 'Majority';
                                case 5
                                    decisiontypestr = 'Median';
                            end
                            vargout.multiclass = [vargout.multiclass ' (' decisiontypestr ')']; 
                        case 2 % ECOC
                            vargout.multiclass = [vargout.multiclass ', Error-correcting output codes'];
                            switch params.MULTI.decoding
                                case 1
                                    decodestr = 'Hamming distance';
                                case 2
                                    decodestr = 'Euclidean distance';
                                case 3
                                    decodestr = 'Laplacian decoding';
                                case 4
                                    decodestr = 'Attenuated euclidean distance';
                                case 5
                                    decodestr = 'Linear loss-based decoding';
                            end
                            vargout.multiclass = [vargout.multiclass ' (' decodestr ')']; 

                        case 3 % Hierarchical One-Vs-One

                        case 4 % DAG


                    end    
                else
                    vargout.multiclass = [vargout.multiclass ': One-vs-One-Max-Wins'];
                end
            else
                vargout.multiclass = 'Disabled.';
            end
        else
            vargout.multiclass = 'NA';
        end

        
    case 'oocv'
        
        if isfield(params.OOCV,'meanflag'), meanflag = params.OOCV.meanflag; end
        if isfield(params.OOCV,'groupmode'), groupmode = params.OOCV.groupmode; end
        if isfield(params.OOCV,'trainwithCV2Ts'), trainwithCV2Ts = params.OOCV.trainwithCV2Ts; end
        if isfield(params.OOCV,'savemodels'), savemodels = params.OOCV.savemodels; end
        if isfield(params.OOCV,'saveoocvdata'), saveoocvdata = params.OOCV.saveoocvdata; end
        
        switch meanflag
            case 1
                vargout.meanflagstr = 'Aggregate all base learnerns into one big ensemble';
            case 2
                vargout.meanflagstr = 'Compute mean decision for each CV1 partition ensemble';
        end

        switch groupmode
            case 1
                vargout.groupmodestr = 'OOCV prediction at binary predictors'' optimum parameters.';
            case 2
                vargout.groupmodestr = 'OOCV prediction at multi-group predictors'' optimum parameters.';
            case 3
                vargout.groupmodestr = 'OOCV prediction at binary & multi-group predictors'' optimum parameters.';

        end

        switch trainwithCV2Ts
            case 1
               vargout.trainwithCV2Tsstr = 'as defined for CV model training';
            case 2
               vargout.trainwithCV2Tsstr = 'use CV1 + CV2 data'; 
        end

end

return

%% IMRELIEF
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function description = Get_IMReliefDescription(params)

if iscell(params.imrelief)
    nclass = numel(params.imrelief);
else
    nclass = 1;
end
description = '';

for curclass = 1:nclass 
    
    if iscell(params.imrelief)
        imrelief = params.imrelief{curclass};
        hdr = sprintf('Cl%g', curclass);
    else
        imrelief = params.imrelief;
        hdr = 'All';
    end
    sigmastr = sprintf('%g ',imrelief.sigma);
    lambdastr = sprintf('%g ',imrelief.lambda);

    description = sprintf('%s%s: %s, s=[%s], l=[%s], i=%g', description, hdr, imrelief.distance, sigmastr, lambdastr, imrelief.maxiter);      
    
    if curclass < nclass, description = sprintf('%s, ', description); end

end

return

%% GFLIP
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function description = Get_GflipDescription(params)
if params.gflip.utilfunc == 3
    if isnumeric(params.gflip.extra_param.beta)
        betastr = sprintf(', beta = %g', ...
            params.gflip.extra_param.beta);
    else
        betastr = ', beta = auto';
    end
else
    betastr = '';
end

description = sprintf('utility = %s, # start points = %g, block size = %g %s', ...
    params.gflip.extra_param.utility, ...
    params.gflip.extra_param.start_points, ...
    params.gflip.extra_param.block_size, ...
    betastr);   
return

%% RGS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function description = Get_RGSDescription(params)
if isfield(params.RGS.extra_param,'beta')
    betastr = num2str(params.RGS.extra_param.beta);
else
    betastr = 'auto';
end
if isfield(params.RGS.extra_param,'k')
    kstr = num2str(params.RGS.extra_param.k);
else
    kstr = 'auto';
end

description = sprintf('# epochs = %g, # starts = %g, # beta = %s, %s neighbors ', ...
    params.RGS.extra_param.epochs, ...
    params.RGS.extra_param.num_starts, ...
    betastr, ...
    kstr);       
return

%% RSS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function description = Get_RSSDescription(params)

description = sprintf('RSS: # nperms = %g, Min. perc. feats = %g, Max. perc. feats = %g', ...
    params.RSS.nperms, ...
    params.RSS.nMinFeatsPerc, ...
    params.RSS.nMaxFeatsPerc);

return

%% FEAST
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function description = Get_FEASTDescription(params, longfl)
if ~exist('longfl','var') || isempty(longfl)
    longfl = 0;
end

if longfl
    switch params.MethodStr
        case 'mim'
            description = 'Mutual Information Maximization (MIM)';
        case 'mrmr'
            description = 'Maximum Relevance Minimum Redundancy (MRMR)';
        case 'cmim'
            description = 'Conditional Mutual Information Maximization (CMIM)';
        case 'jmi'
            description = 'Joint Mutual Information (JMI)';
        case 'disr'
            description = 'Double Input Symmetrical Relevance (DISR)';
        case 'cife'
            description = 'Conditional Infomax Feature Extraction (CIFE)';
        case 'icap'
            description = 'Interaction Capping (ICAP)';
        case 'condred'
            description = 'CONDRED';
        case 'cmi'
            description = 'CMI';
        case 'mifs'
            description = 'Mutual Information Feature Selection (MIFS)';
    end
else
    description = sprintf('%s', params.FEAST.MethodStr);
end
return
