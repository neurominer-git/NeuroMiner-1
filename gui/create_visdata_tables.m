function VTBL = create_visdata_tables(v, VTBL, ind, act)

if ~exist('act','var') || isempty(act)
    error('Action string is required')
elseif strcmp(act,'reorder') && ( (~exist('VTBL','var') || isempty(VTBL)) || (~exist('ind','var') || isempty(ind)) )
    error('Table data or sorting vector missing')
end

switch act
    case 'create'
        if isempty(VTBL)
            VTBL.tbl = struct('array',[], 'colnames', [], 'rownames',[], 'ind', []);
        end
        VTBL.params = v.params;
        if iscell(v.MEAN)
            if isempty(ind), ind = (1:size(v.MEAN{1},1))'; end
            nclass = numel(v.MEAN);
            for curclass = 1:nclass
                VTBL.tbl(curclass).ind = ind;
                VTBL.tbl(curclass).array = [ind v.MEAN{curclass}(ind), v.SE{curclass}(ind), ...
                    v.MEAN_CV2{curclass}(ind), v.SE_CV2{curclass}(ind), ...
                    v.CVRatio{curclass}(ind), v.CVRatio_CV2{curclass}(ind), v.Prob_CV2{curclass}(ind)] ;
                VTBL.tbl(curclass).colnames = ...
                    {'Feature', 'SortIndex', 'Mean_W', 'StErr_W', 'Mean_W_GM', 'StErr_W_GM', 'CVratio', 'CVratio_GM', 'P_Relia95CI_GM'};
                if isfield(v,'SignBased_CV2')
                    if iscell(v.SignBased_CV2)
                        vx = v.SignBased_CV2{curclass}(ind);
                    else
                        vx = v.SignBased_CV2(ind);
                    end
                    VTBL.tbl(curclass).array = [ VTBL.tbl(curclass).array vx, v.SignBased_CV2_z{curclass}(ind), v.SignBased_CV2_p_uncorr{curclass}(ind), v.SignBased_CV2_p_fdr{curclass}(ind) ];
                    VTBL.tbl(curclass).colnames = [ VTBL.tbl(curclass).colnames 'SignConst', 'Z_SignConst', 'P_SignConst', 'Pfdr_SignConst' ];
                end
                VTBL.tbl(curclass).rownames = v.params.features(ind)';
                if isfield(v,'Spearman_CV2')
                    VTBL.tbl(curclass).array = [VTBL.tbl(curclass).array  v.Spearman_CV2{curclass}(ind)];
                    VTBL.tbl(curclass).colnames = [VTBL.tbl(curclass).colnames, 'Spearman_GM'];     
                    VTBL.tbl(curclass).array = [VTBL.tbl(curclass).array  v.Pearson_CV2{curclass}(ind)];
                    VTBL.tbl(curclass).colnames = [VTBL.tbl(curclass).colnames, 'Pearson_GM'];
                    if isfield(v,'Spearman_CV2_p_uncorr')
                        VTBL.tbl(curclass).array = [VTBL.tbl(curclass).array  v.Spearman_CV2_p_uncorr{curclass}(ind)];
                        VTBL.tbl(curclass).colnames = [VTBL.tbl(curclass).colnames, 'P_Spearman_GM'];     
                        VTBL.tbl(curclass).array = [VTBL.tbl(curclass).array  v.Pearson_CV2_p_uncorr{curclass}(ind)];
                        VTBL.tbl(curclass).colnames = [VTBL.tbl(curclass).colnames, 'P_Pearson_GM'];
                        VTBL.tbl(curclass).array = [VTBL.tbl(curclass).array  v.Spearman_CV2_p_fdr{curclass}(ind)];
                        VTBL.tbl(curclass).colnames = [VTBL.tbl(curclass).colnames, 'Pfdr_Spearman_GM'];     
                        VTBL.tbl(curclass).array = [VTBL.tbl(curclass).array  v.Pearson_CV2_p_fdr{curclass}(ind)];
                        VTBL.tbl(curclass).colnames = [VTBL.tbl(curclass).colnames, 'Pfdr_Pearson_GM'];
                    end
                end
                if isfield(v, 'FeatProb')
                    VTBL.tbl(curclass).array = [VTBL.tbl(curclass).array  v.FeatProb{1}(ind,curclass)];
                    VTBL.tbl(curclass).colnames = [VTBL.tbl(curclass).colnames, 'FeatSelProb'];
                end
                if isfield(v, 'PermZ_CV2')
                    VTBL.tbl(curclass).array = [VTBL.tbl(curclass).array  v.PermZ_CV2{curclass}(ind)];
                    VTBL.tbl(curclass).colnames = [VTBL.tbl(curclass).colnames, 'Z_Perm'];
                    VTBL.tbl(curclass).array = [VTBL.tbl(curclass).array  v.PermProb_CV2{curclass}(ind)];
                    VTBL.tbl(curclass).colnames = [VTBL.tbl(curclass).colnames, 'P_PermProb'];
                    VTBL.tbl(curclass).array = [VTBL.tbl(curclass).array  v.PermProb_CV2_FDR{curclass}(ind)];
                    VTBL.tbl(curclass).colnames = [VTBL.tbl(curclass).colnames, 'Pfdr_PermProb'];
                end
            end
        else
            if isempty(ind), ind = (1:size(v.MEAN,1))'; end
            VTBL.tbl.ind = ind; 
            VTBL.tbl.array = [ind v.MEAN(ind), v.SE(ind), v.MEAN_CV2(ind), v.SE_CV2(ind), v.CVRatio(ind), v.CVRatio_CV2(ind), v.Prob_CV2(ind), ...
                v.SignBased_CV2(ind), v.SignBased_CV2_z(ind), v.SignBased_CV2_p_uncorr(ind), v.SignBased_CV2_p_fdr(ind) ];
            VTBL.tbl.colnames = {'Feature', 'SortIndex', 'Mean_W', 'StErr_W', 'Mean_W_GM', 'StErr_W_GM', 'CVratio', 'CVratio_GM', 'P_Relia95CI_GM', ...
                'SignConst', 'Z_SignConst', 'P_SignConst', 'Pfdr_SignConst'};
            VTBL.tbl.rownames = v.params.features(ind)';
            if isfield(v,'Spearman_CV2') && isfield(v,'Spearman_CV2_uncorr')
                VTBL.tbl.array = [VTBL.tbl.array  v.Spearman_CV2(ind), v.Pearson_CV2(ind), v.Spearman_CV2_p_uncorr(ind), v.Pearson_CV2_p_uncorr(ind), v.Spearman_CV2_p_fdr(ind), v.Pearson_CV2_p_fdr(ind)];
                VTBL.tbl.colnames = [VTBL.tbl.colnames, 'Spearman_GM', 'Pearson_GM', 'P_Spearman_GM', 'P_Pearson_GM','Pfdr_Spearman_GM', 'Pfdr_Pearson_GM'];           
            end
            if isfield(v, 'FeatProb')
                VTBL.tbl.array = [VTBL.tbl.array  v.FeatProb(ind)];
                VTBL.tbl.colnames = [VTBL.tbl.colnames, 'FeatSelProb'];
            end
            if isfield(v, 'PermZ_CV2')
                VTBL.tbl.array = [VTBL.tbl.array  v.PermZ_CV2(ind) v.PermProb_CV2(ind) v.PermProb_CV2_FDR(ind)];
                VTBL.tbl.colnames = [VTBL.tbl.colnames, 'Z_Perm', 'P_Perm_CV2', 'Pfdr_Perm_CV2'];
            end
        end
        
    case 'reorder'
        
        if numel(VTBL.tbl)>1
            nclass = numel(VTBL.tbl);
            for curclass=1:nclass
                VTBL.tbl(curclass).ind = ind;
                VTBL.tbl(curclass).array = VTBL.tbl(curclass).array(ind,:);
                VTBL.tbl(curclass).rownames = VTBL.tbl(curclass).rownames(ind);
            end
        else
            VTBL.tbl.ind = ind;
            VTBL.tbl.array = VTBL.tbl.array(ind,:);
            VTBL.tbl.rownames = VTBL.tbl.rownames(ind);
        end
end