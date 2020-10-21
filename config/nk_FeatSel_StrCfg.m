function [strtype, suff] = nk_FeatSel_StrCfg(FEATSEL, varind)

suffx='';
% Build output data strings
if isfield(FEATSEL,'ACTPARAM')
    
    for i=1:numel(FEATSEL.ACTPARAM)
        
        if ~isempty(FEATSEL.ACTPARAM{i})
            
            switch FEATSEL.ACTPARAM{i}.cmd
                
                case 'scale'
                    suffx = [ suffx '_scl' ];
                case 'standardize'
                    suffx = [ suffx '_std' ];
                case 'discretize'
                    suffx = [ suffx '_dsc' ];
                case 'symbolize'
                    suffx = [ suffx '_sym' ];
                case 'correctnuis'
                    suffx = [ suffx '_nuis' ];
                    for j=1:numel(FEATSEL.ACTPARAM{i}.COVAR)
                        suffx = [ suffx '-' num2str(FEATSEL.ACTPARAM{i}.COVAR(j)) ];
                    end
                case 'norm'
                    suffx = [ suffx '_norm' ];
                    for j=1:numel(PREPROC.ACTPARAM{i}.IND)
                        suffx = [ suffx '-' num2str(PREPROC.ACTPARAM{i}.IND(j)) ];
                    end
                case 'rankfeat'
                    switch FEATSEL.ACTPARAM{i}.weightmethod
                        case 1
                            suffx = [ suffx '_upw-' FEATSEL.ACTPARAM{i}.labeldesc '-' FEATSEL.ACTPARAM{i}.algostr];
                        case 2
                            suffx = [ suffx '_dnw-' FEATSEL.ACTPARAM{i}.labeldesc '-' FEATSEL.ACTPARAM{i}.algostr];
                    end
            end
        end
    end
end

switch FEATSEL.varflag
    case 1
        switch FEATSEL.ranktype 
            case 1, 
                strtype = ['pearson' suffx];
            case 2
                strtype = ['FDR' suffx];
            case 3
                strtype = ['FTest' suffx];
            case 4
                strtype = 'MI';               
        end
    case 2
        switch FEATSEL.ranktype 
            case 1, 
                strtype = ['PLS' suffx];
            case 2
                strtype = ['kPLS' suffx];
            case 3
                strtype = ['PCA' suffx];
        end
    case 3
        switch FEATSEL.ranktype
            case 1
                strtype = ['ams' suffx];
            case 2
                switch FEATSEL.simba.utilfunc
                    case 1
                        strtype = ['simba-lin' suffx];
                    case 2
                        strtype = ['simba-sigm' suffx];
                end
            case 3
                switch FEATSEL.gflip.utilfunc
                    case 1
                        strtype = ['gflip-zo' suffx];
                    case 2
                        strtype = ['gflip-lin' suffx];
                    case 3
                        strtype = ['gflip-sigm' suffx];
                end
            case 4
                switch FEATSEL.gflip.utilfunc
                    case 1
                        strtype = ['gflip-zo' suffx];
                    case 2
                        strtype = ['gflip-lin' suffx];
                    case 3
                        strtype = ['gflip-sigm' suffx];
                end
                switch FEATSEL.simba.utilfunc
                    case 1
                        strtype = ['simba-lin_' strtype];
                    case 2
                        strtype = ['simba-sigm_' strtype];
                end
            case 5
                strtype = ['RGS' suffx];
            case 6
                strtype = ['IMRelief' suffx];
            case 7
                strtype = [FEATSEL.FEAST.MethodStr suffx];
            case 8
                strtype = ['SVM' suffx];
            case 9
                strtype = ['ccSVM' suffx];
            case 10
                strtype = ['AUC' suffx ];
            case 11
                strtype = ['Pearson' suffx ];
            case 12
                strtype = ['Spearman' suffx ];
            case 13
                strtype = ['FScore' suffx ];
            case 14
                strtype = ['Relief' suffx ];
        end
end

switch FEATSEL.cubetype
    case 1
        strtype = [strtype '_sp4x4'];
    case 2
        strtype = [strtype '_sp9x3'];
    case 3
        strtype = [strtype '_sp-fwhm' num2str(FEATSEL.cubefwhm)];
end

switch FEATSEL.varflag
    case 1
        if ~FEATSEL.permflag 
            permsuff = '';
        else
            permsuff = ['_p' num2str(FEATSEL.nperms)];
        end
        switch FEATSEL.cv,
            case 0
                suff = permsuff;
            case 1
                suff = ['_loo' permsuff];
            case 2
                suff = ['_cv' num2str(FEATSEL.kfold) permsuff];
            case 3
                suff = ['_bs' num2str(FEATSEL.nboot) permsuff ];
        end
    case 2
        suff = [];
        if FEATSEL.bootflag 
            suff = ['_bs' num2str(FEATSEL.nboot)];
        end
        if FEATSEL.permflag
            suff = [suff '_p' num2str(FEATSEL.nperms)];
        end
        
    case 3
        suff='';
        
end

for i=1:numel(varind)
    suff = [suff '_var' num2str(varind(i))];
end

return