function [I2, I1, filefound] = nk_VisXHelper(act, nM, nclass, decompfl, permfl, ix, jx, I2, inp, ll, nperms, I1)
global FUSION SVM

linsvmfl = determine_linsvm_flag(SVM);
filefound = false;
switch act
    
    case 'init'
        % Initialize arrays in data containers for CV2-level statistics
        if any(permfl)
            I2.VCV2PERM      = cell(nclass,nM);                  % Container for feature space null distribution estimations: p value (each cell for one modality)
            I2.VCV2PERM_FDR  = cell(nclass,nM);                  % Container for feature space null distribution estimations: p value (each cell for one modality)
            I2.VCV2ZSCORE    = cell(nclass,nM);                  % Container for feature space null distribution estimations: Z score (each cell for one modality)
            I2.VCV2MPERM_CV2 = nan(nclass,ix*jx);
            I2.VCV2MPERM     = nan(nclass,ix*jx);
        end
        if any(~decompfl)
            I2.VCV2PEARSON              = cell(nclass,nM);                  % Container for feature spaces' univariate pearson correlation coefficients
            I2.VCV2SPEARMAN             = cell(nclass,nM);                  % Container for feature spaces' univariate spearman correlation coefficients
            I2.VCV2PEARSON_UNCORR_PVAL  = cell(nclass,nM);
            I2.VCV2SPEARMAN_UNCORR_PVAL = cell(nclass,nM);
            I2.VCV2PEARSON_FDR_PVAL     = cell(nclass,nM);
            I2.VCV2SPEARMAN_FDR_PVAL    = cell(nclass,nM);
            if linsvmfl, 
                I2.VCV2PVAL_ANALYTICAL = cell(nclass,nM);
                I2.VCV2PVAL_ANALYTICAL_FDR = cell(nclass,nM);
            end
        end
        I2.VCV2VCV1         = cell(nclass,nM);                  % Container for concatenation of CV1 feature spaces' 
        I2.VCV2SUM          = cell(nclass,nM);                  % Container for feature spaces' sum values of CV-generated feature weights
        I2.VCV2SQ           = cell(nclass,nM);                  % Container for feature spaces' sqrt values of CV-generated feature weights
        I2.GCV2SUM          = cell(nclass,nM);
        I2.VCV2MEAN         = cell(nclass,nM);                  % Container for feature spaces' mean vectors of CV-generated feature weights
        I2.VCV2STD          = cell(nclass,nM);                  % Container for feature spaces' std vectors of CV-generated feature weights
        I2.VCV2PROB         = cell(nclass,nM);                  % Container for feature spaces' reliability masks of CV-generated feature weights (reliable means |mean| > sterr)
        I2.VCV2SEL          = cell(nclass,nM);                  % Vectors with #number of models per feature info
        I2.PCV2SUM          = cell(nclass,nM);                  % Vectors with 'times selected' feature info
        I2.VCV2NMODEL       = zeros(nclass,1);
    
    case 'initI1'
        I1.VCV1NMODEL       = zeros(nclass,1);
        I1.VCV1SUM          = cell(nclass,nM);
        I1.VCV1SQ           = cell(nclass,nM);
        I1.VCV1MEAN         = cell(nclass,nM);
        I1.VCV1STD          = cell(nclass,nM);
        if any(~decompfl)
            I1.PCV1SUM                  = cell(nclass,nM);
            I1.VCV1PEARSON              = cell(nclass,nM);
            I1.VCV1SPEARMAN             = cell(nclass,nM);
            I1.VCV1PEARSON_UNCORR_PVAL  = cell(nclass,nM);
            I1.VCV1SPEARMAN_UNCORR_PVAL = cell(nclass,nM);
            I1.VCV1PEARSON_FDR_PVAL     = cell(nclass,nM);
            I1.VCV1SPEARMAN_FDR_PVAL    = cell(nclass,nM);
            if linsvmfl, 
                I1.VCV1PVAL_ANALYTICAL = cell(nclass,nM);
                I1.VCV1PVAL_ANALYTICAL_FDR = cell(nclass,nM);
            end
        end
        if any(permfl)
            I1.VCV1PERM     = cell(nclass,nM);
            I1.VCV1PERM_FDR = cell(nclass,nM);
            I1.VCV1ZSCORE   = cell(nclass,nM);
            I1.VCV1WPERM    = cell(nclass,1);
            I1.VCV1MPERM    = cell(nclass,1);
        end
        I2 = [];
        
    case 'accum'
        
        if exist('I1','var')
            if ischar(I1) && exist(I1,'file'), 
                try
                    load(I1); 
                    filefound = true;
                catch
                    cprintf('red','\nCould not load file %s ', I1);
                    filefound = false;
                    return;
                end
            end
        else
            error('Visualization structure I1 has to be provided!')
        end

        for h=1:nclass

            if any(permfl), 
                I2.VCV2MPERM(h,ll) 		  = mean(I1.VCV1MPERM{h});
                I2.VCV2MPERM_CV2(h,ll) 	  = sum(I1.VCV1MPERM_CV2{h})/nperms(1);
            end

            % Loop through modalities
            for n = 1:nM

                [D,~,~,badcoords]         = getD(FUSION.flag, inp, n); badcoords = ~badcoords; 
                if ~decompfl(n),			
                    if isempty(I2.PCV2SUM{h,n}), 
                        I2.PCV2SUM{h, n}(badcoords)        = I1.PCV1SUM{h, n}(badcoords);
                    else
                        I2.PCV2SUM{h, n}(badcoords)        = I2.PCV2SUM{h, n}(badcoords)' + I1.PCV1SUM{h, n}(badcoords);
                    end
                    % Changed from median to mean computation (30.12.2018)
                    % log-transform first then compute the mean for P values 
                    I2.VCV2PEARSON{h, n}                = [ I2.VCV2PEARSON{h, n}                nm_nanmean(I1.VCV1PEARSON{h, n},2) ];
                    I2.VCV2SPEARMAN{h, n}               = [ I2.VCV2SPEARMAN{h, n}               nm_nanmean(I1.VCV1SPEARMAN{h, n},2) ];
                    I2.VCV2PEARSON_UNCORR_PVAL{h, n}    = [ I2.VCV2PEARSON_UNCORR_PVAL{h, n}    nm_nanmean(-log10(I1.VCV1PEARSON_UNCORR_PVAL{h, n}),2) ];
                    I2.VCV2SPEARMAN_UNCORR_PVAL{h, n}   = [ I2.VCV2SPEARMAN_UNCORR_PVAL{h, n}   nm_nanmean(-log10(I1.VCV1SPEARMAN_UNCORR_PVAL{h, n}),2) ];
                    I2.VCV2PEARSON_FDR_PVAL{h, n}       = [ I2.VCV2PEARSON_FDR_PVAL{h, n}       nm_nanmean(-log10(I1.VCV1PEARSON_FDR_PVAL{h, n}),2) ];
                    I2.VCV2SPEARMAN_FDR_PVAL{h, n}      = [ I2.VCV2SPEARMAN_FDR_PVAL{h, n}      nm_nanmean(-log10(I1.VCV1SPEARMAN_FDR_PVAL{h, n}),2) ];
                    if linsvmfl && isfield(I1,'VCV1PVAL_ANALYTICAL')
                        I2.VCV2PVAL_ANALYTICAL{h, n}    = [ I2.VCV2PVAL_ANALYTICAL{h, n}        nm_nanmean(-log10(I1.VCV1PVAL_ANALYTICAL{h,n}),2) ];
                        I2.VCV2PVAL_ANALYTICAL_FDR{h, n}= [ I2.VCV2PVAL_ANALYTICAL_FDR{h, n}    nm_nanmean(-log10(I1.VCV1PVAL_ANALYTICAL_FDR{h,n}),2) ];
                    end
                end
                if any(permfl)
                    % Changed from median to mean computation for P values (30.12.2018)
                    if isempty(I2.VCV2PERM{h, n})
                        I2.VCV2PERM{h, n}       = double(nm_nanmean(I1.VCV1PERM{h, n},2)) ;
                        I2.VCV2PERM_FDR{h, n}   = double(nm_nanmean(I1.VCV1PERM_FDR{h, n},2));
                        I2.VCV2ZSCORE{h, n}     = double(nm_nansum(I1.VCV1ZSCORE{h, n},2));
                    else
                        I2.VCV2PERM{h, n}       = [ I2.VCV2PERM{h,n} double(nm_nanmean(I1.VCV1PERM{h, n},2)) ];
                        I2.VCV2PERM_FDR{h, n}   = [ I2.VCV2PERM_FDR{h,n} double(nm_nanmean(I1.VCV1PERM_FDR{h, n},2)) ];
                        I2.VCV2ZSCORE{h, n}     = [ I2.VCV2ZSCORE{h,n} double(nm_nansum(I1.VCV1ZSCORE{h, n},2)) ];
                    end
                end 
                I1.VCV1NMODEL(h)                = size(I1.VCV1{h,n},2);
                I1.VCV1SEL{h,n}                 = sum(~isnan(I1.VCV1{h,n}),2);
                I1.VCV1MEAN{h,n}                = nm_nanmedian(I1.VCV1{h,n},2);
                I1.VCV1SUM{h,n}                 = nm_nansum(I1.VCV1{h,n},2);
                I1.VCV1SQ{h,n}                  = nm_nansum(I1.VCV1{h,n}.^2,2);
                I1.VCV1STD{h,n}                 = (nm_nanstd(I1.VCV1{h,n},2)./sqrt(I1.VCV1SEL{h, n}))*1.96;
                indMEANgrSE                     = abs(I1.VCV1MEAN{h,n}) > I1.VCV1STD{h,n};
                if isempty(I2.VCV2SUM{h, n})
                    I2.GCV2SUM{h, n}            = nan(D,ix*jx,'double');
                    I2.VCV2PROB{h, n}           = indMEANgrSE;
                    I2.VCV2SUM{h, n}            = I1.VCV1SUM{h,n};
                    I2.VCV2SUM2{h,n}            = (I1.VCV1SUM{h,n}.^2)./I1.VCV1SEL{h, n};
                    I2.VCV2SQ{h, n}             = I1.VCV1SQ{h,n};
                    I2.VCV2MEAN{h, n}           = I1.VCV1MEAN{h,n};
                    I2.VCV2STD{h, n}            = I1.VCV1STD{h,n};
                    I2.VCV2SEL{h, n}            = I1.VCV1SEL{h, n};
                    I2.VCV2VCV1{h,n}            = I1.VCV1{h,n};
                else
                    I2.GCV2SUM{h, n}(:,ll)      = I1.VCV1SUM{h, n} ./  I1.VCV1SEL{h,n};  
                    I2.VCV2PROB{h, n}           = [I2.VCV2PROB{h, n}   indMEANgrSE ];
                    I2.VCV2SUM{h, n}            = [I2.VCV2SUM{h, n}    I1.VCV1SUM{h, n}];
                    I2.VCV2SUM2{h,n}            = [I2.VCV2SUM2{h, n}   (I1.VCV1SUM{h,n}.^2)./I1.VCV1SEL{h, n}];
                    I2.VCV2SQ{h, n}             = [I2.VCV2SQ{h, n}     I1.VCV1SQ{h, n} ];
                    I2.VCV2MEAN{h, n}           = [I2.VCV2MEAN{h, n}   I1.VCV1MEAN{h, n}];
                    I2.VCV2STD{h, n}            = [I2.VCV2STD{h, n}    I1.VCV1STD{h, n}];
                    I2.VCV2SEL{h, n}            = I2.VCV2SEL{h, n}     + I1.VCV1SEL{h, n};
                    I2.VCV2VCV1{h,n}            = [I2.VCV2VCV1{h,n}    I1.VCV1{h,n}];
                end 
            end
            I2.VCV2NMODEL(h) = I2.VCV2NMODEL(h) + I1.VCV1NMODEL(h);
        end

end
