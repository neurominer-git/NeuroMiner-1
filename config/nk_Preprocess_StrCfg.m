function strout = nk_Preprocess_StrCfg(ix, jx)

global PREPROC

if isfield(PREPROC,'FEATSEL') && ...
        isfield(PREPROC.FEATSEL,'BINMOD') && ...
        PREPROC.FEATSEL.BINMOD ==1 && ...
        PREPROC.BINMOD == 0;
    suffx = '_Concat';
else
    suffx='';
end

if isfield(PREPROC,'SPATIAL')
    switch PREPROC.SPATIAL.cubetype
        case 4
            suffx = [suffx '_flt-FWHM'];
        case 3
            suffx = [suffx '_flt-v27n'];
        case 2
            suffx = [suffx '_flt-r4n'];
    end
end

if isfield(PREPROC,'ACTPARAM')
    
    for i=1:numel(PREPROC.ACTPARAM)
        
        if ~isempty(PREPROC.ACTPARAM{i})
            
            switch PREPROC.ACTPARAM{i}.cmd
                case 'selectfeat'
                    if ~isfield(PREPROC.ACTPARAM{i},'MODE')
                        MODE = 1;
                    else
                        MODE = PREPROC.ACTPARAM{i}.MODE;
                    end
                    switch MODE
                        case 1
                            if PREPROC.ACTPARAM{i}.CLUST == 1
                                suffx = [suffx '_cl'];
                            else
                                suffx = [suffx '_vx'];
                            end
                        case 2
                            suffx = [suffx '_wg'];
                    end
                case 'scale'
                    suffx = [ suffx '_scl' ];
                case 'standardize'
                    switch PREPROC.ACTPARAM{i}.METHOD
                        case {'standardization', 'standardization using median'}
                            stdsuff = '';
                        case 'standardization using mean'
                            stdsuff = '-mn';
                        case 'mean-centering' 
                            stdsuff = '-mc';
                        case 'l1-median centering'
                            stdsuff = '-l1-mc';
                        case 'qn-standardization'
                            stdsuff = '-qn';
                        case 'sn-standardization'
                            stdsuff = '-sn';
                    end
                    suffx = [ suffx '_std' stdsuff  ];
                case 'discretize'
                    suffx = [ suffx '_dsc' ];
                case 'symbolize'
                    suffx = [ suffx '_sym' ];
                case 'correctnuis'
                    suffx = [ suffx '_nuis' ];
                    for j=1:numel(PREPROC.ACTPARAM{i}.COVAR)
                        suffx = [ suffx '-' num2str(PREPROC.ACTPARAM{i}.COVAR(j)) ];
                    end
                case 'normalize'
                    suffx = [ suffx '_norm' ];
                    for j=1:numel(PREPROC.ACTPARAM{i}.IND)
                        suffx = [ suffx '-' num2str(PREPROC.ACTPARAM{i}.IND(j)) ];
                    end
                case 'reducedim'
                    suffx = [suffx '_' PREPROC.ACTPARAM{i}.DR.RedMode];
                case 'rankfeat'
                    if sum(any(PREPROC.ACTPARAM{i+1}.W_ACT.threshvec))>0
                        if PREPROC.ACTPARAM{i+1}.W_ACT.clustflag
                            suffx = [suffx '_rnk-' PREPROC.ACTPARAM{i}.RANK.algostr '-hcl'];
                        else
                            suffx = [suffx '_rnk-' PREPROC.ACTPARAM{i}.RANK.algostr '-hvx'];
                        end
                    else
                        suffx = [suffx '_rnk-' PREPROC.ACTPARAM{i}.RANK.algostr '-swg'];
                    end
                case 'devmap'
                     suffx = [suffx '-devmp-' PREPROC.ACTPARAM{i}.DEVMAP.algostr ];
            end
        end
    end
end

if (exist('ix','var') && exist('jx','var')) && (~isempty(ix) && ~isempty(jx)) 
    strout = [suffx '_oCV' num2str(ix) '.' num2str(jx)]; 
else
    strout = suffx;
end

return