function T = printvis_ranked_cross_binary(analysis, fld, reorder, features, nfeat, act)
% T = printvis_ranked_cross_binary(analysis, fld, reorder, features, nfeat, act)
% ==============================================================================
% The algorithm computes the mean of fld across available binary
% classifiers, ranks amd prints ranked features according to mean(fld).
%
% Input arguments:
% analysis:     NM.analysis{x} structure
% fld           can be one of the following depending on your visualisation 
%               results. Check NM.analysis{x}.visdata to now exactly what 
%               is available
%   MEAN_CV2
%   SE_CV2: 
%   CVRatio_CV2: 
%   Prob_CV2: 
%   Pearson_CV2: 
%   Spearman_CV2: 
%   Pearson_CV2_p_uncorr: 
%   Spearman_CV2_p_uncorr: 
%   Pearson_CV2_p_fdr: 
%   Spearman_CV2_p_fdr: 
%   Pearson_CV2_STD: 
%   Spearman_CV2_STD: 
%   Pearson_CV2_p_uncorr_STD: 
%   Spearman_CV2_p_uncorr_STD: 
%   Analytical_p: 
%   Analyitcal_p_fdr: 
%   PermProb_CV2: 
%   PermProb_CV2_FDR: 
%   PermZ_CV2:
%
% reorder:      reorder features in plot (set to [] for default operation)
% features:     cell array of strings with names of your features
% nfeat:        top nfeat features to be printed
% act:          'print' or 'info' (use the latter to obtain info on what is
%               available in the visualization structure)
%
% Output arguments:
% T             Table structure containing the visualized data
% =========================================================================
% 12/2018 Nikolaos Koutsouleris

nclass = numel(analysis.visdata{1}{1}.MEAN);

if ~exist('act','var') || isempty(act), act = 'info'; end
switch act
    case 'print'
        if ~exist('fld','var') || isempty(fld), fld = 'MEAN_CV2'; end
        if ~exist('reorder','var'), reorder = []; end
        V = [];
        try
            for i=1:nclass
                if size(analysis.visdata{1}{1}.(fld),1)>1
                    V = [V analysis.visdata{1}{1}.(fld){i}];
                else
                    V = [V analysis.visdata{1}{1}.(fld){1}(:,i)];
                end
            end
        catch
            cprintf('red','Requested visualization parameter not found.')
            T = printvis_ranked_cross_binary(analysis, [], [], [], [], 'info');
            return
        end
        if ~exist('nfeat', 'var') || isempty(nfeat), nfeat = size(V,1); end 
        if isempty(reorder)
            mV = mean(V,2);
            [~,ind] = sort(mV, 'descend');
            V = V(ind,:); features = features(ind);
        end
        V = V(1:nfeat,:); features = features(1:nfeat);
        figure;ha=tight_subplot(1,nclass);
        Ylimits = [0.5 nfeat+0.5];
        if min(V(:))>0,
            Xlimits = [0 max(V(:))];
        else
            Xlimits = [min(V(:)) max(V(:))];
        end
        for i=1:nclass
            tit = analysis.params.cv.class{1,1}{i}.groupdesc;
            barh(ha(i),V(:,i));
            ha(i).YLim = Ylimits;
            ha(i).XLim = Xlimits;
            ha(i).Title.String = tit;
            if i==1, 
                ha(1).YTickLabel = features; 
            else
                ha(i).YTickLabel = '';
            end
        end
        T = [cell2table(features') array2table(V)];
    case 'info'
        fprintf('\n%g Classifiers detected in analysis ''%s''.',nclass, analysis.id);
        fprintf('\nThe analysis structure contains the following visualization fields:')
        flds = fieldnames(analysis.visdata{1}{1});
        for i=1:numel(flds)
            fprintf('\n\t%s',flds{i});
        end
        fprintf('\n');
        T=flds;
end
