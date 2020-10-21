%==========================================================================
%FORMAT AgreeFeat = ProbabilisticFea( F, EnsStrat )                                    
%==========================================================================
%Determines optimum feature space through the between-classifier agreement
%of selected features.
%
%This function works in three different operation modes:
%-------------------------------------------------------
%(1) Features are select based on a threshold of between-CV1-feature mask
%agreements. If the user-defined threshold does not return any feature
%meeting the minimum number of features, a user-defined tolerance window 
%is systematically tested. An error is generated if this testing does 
%result in an empty feature mask.
%(2) The features are sorted according to their cross-CV1 agreement and
%the top N features are returned.
%(3) The features are sorted according to their cross-CV1 agreement and
%the top %N features are returned.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(c) Nikolaos Koutsouleris, 09 / 2015
function [AgreeFeat, AgreePerc] = ProbabilisticFea(F, EnsStrat, q)
global VERBOSE
%persistent h_probabilisticfea

% Compute the features' between-classifier conistency
% as percentage of between-classifier "agreement"
[nF, mF] = size(F);
if ~exist('q','var') || isempty(q), q=1; end
AgreePerc = sum(F~=0,2)*100./mF;

switch EnsStrat.Mode
    case 1
        flg = false;
        for i = EnsStrat.Perc(q) : -1 : EnsStrat.Perc(q) - EnsStrat.TolWin
            tF = AgreePerc >= i;
            if sum(tF) >= EnsStrat.MinNum, flg = true; break; end
        end
        if ~flg
            error(['\nNM did not find k = %g suprathreshold features identified within the specified threshold of T = %g (-%g).'...
                    '\nCheck the settings of your probabilistic feature extraction setup!'], ...
                EnsStrat.MinNum, EnsStrat.Perc(q), EnsStrat.TolWin); 
        end
    case 2
        % Absolute
        [~,ind] = sort(AgreePerc,'descend');
        tF = false(nF,1); tF(ind(1:EnsStrat.Perc(q))) = true;
    case 3
        % Percentage selected
        [~,ind] = sort(AgreePerc,'descend');
        thresh = ceil(sum(AgreePerc>0)/100*(100-EnsStrat.Perc(q)));
        tF = false(nF,1); tF(ind(1:thresh)) = true;
        %titlefig = sprintf('Percentage selected: Threshold = %g%%', thresh);
end

if isfield(EnsStrat,'PruneFlag')
    switch EnsStrat.PruneFlag
        case 1
            % Here we prune unselected features for each individual feature
            % mask
            AgreeFeat = logical(F);
            AgreeFeat(bsxfun(@ne,AgreeFeat, tF)) = false;
        case 2
            % Here we take over only the selected features into a single
            % mask for all CV1 training partitions
             AgreeFeat = tF;
    end
else
    AgreeFeat = tF; 
end

if VERBOSE
    h_temp = findobj('Tag','FeatureMatrix');
    if isempty(h_temp)
        h_temp = figure('Name', 'Probabilistic Feature Exraction', ...
            'NumberTitle','off', ...
            'Tag','FeatureMatrix', ...
            'MenuBar','none', ...
            'Color', [0.9 0.9 0.9]);
        h_ax1 = axes(h_temp,'Tag','Orig'); h_ax1.Position(1) = 0.1; h_ax1.Position(3) = 0.25;
        h_ax2 = axes(h_temp,'Tag','Agreement'); h_ax2.Position(1) = 0.4; h_ax2.Position(3) = 0.25; 
        h_ax3 = axes(h_temp,'Tag','Selected'); h_ax3.Position(1) = 0.7; h_ax3.Position(3) = 0.25; 
    else
        h_ax1 = h_temp.Children(1);
        h_ax2 = h_temp.Children(2);
        h_ax3 = h_temp.Children(3);
    end
    imagesc(h_ax1, ~logical(F)); colormap(h_ax1,'gray'); title(h_ax1,'Original Feature matrix'); 
    barh(h_ax2, AgreePerc); set(h_ax2, 'Ydir', 'reverse'); ylim(h_ax2,[0.5 nF-0.5]); title(h_ax2,'Cross CV1 Feature Agreement');
    imagesc(h_ax3,~AgreeFeat); colormap(h_ax3,'gray'); title(h_ax3,'Selected Features'); 
    drawnow
end

