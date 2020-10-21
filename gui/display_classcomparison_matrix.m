function Mx = display_classcomparison_matrix(pvals, ticklabels, mw_all, sw_all, T, G, hlinepos, perftitle)
% =========================================================================
% function Mx = display_classcomparison_matrix(pvals, ticklabels, mw, sw)
% =========================================================================
% Displays P value matrix of pair-wise classifier comparisons 
% and mean(sd) values of classifiers' performance
% 
% Inputs:
% -------
% pvals :           P value matrix ( Pairwise n_class x n_class array )
% ticklabels :      Classifier descriptions (1 x n_class cell arry)
% mw_all :          mean CV2 performance of classifiers (1 x n_class or 1 x n_class x n_comparisons matrix )
% sw_all :          sd CV2 performance of classifiers (1 x n_class, 2 x n_class or 2 x n_class x n_comparisons cube)
% T :               C value scaling for P value matrix
% G :               Performance scores( n_pairs x n_class matrix )
% hlinepos :        optional horizontal line 
% perftitle :       optional performance label

% Outputs:
% --------
% Mx :              discretized classifier matrix
% _________________________________________________________________________
% (c) Nikolaos Koutsouleris, 06/2020

sz = get(0,'ScreenSize');
win_wdth = sz(3)/2; win_hght = sz(4)*0.8; win_x = sz(3)/2 - win_wdth/2; win_y = sz(4)/2 - win_hght/2;
if ~exist('T','var') || isempty(T)
    T = -log10([0.05 0.01 0.001 0.0001 0.00001 0.000001 0.0000001 realmin]) ;
end

if istable(pvals)
    pval_mat_all = pvals.Variables;
else
    pval_mat_all = pvals;
end
if ~exist('perftitle','var') || isempty(perftitle),
    perftitle = 'Performance';
end
nT = numel(T);    
nM = size(pval_mat_all,3);
figure('Name','Comparative Model Analysis', 'Position',[win_x win_y win_wdth win_hght]);         
ax = tight_subplot(2,nM,0.03,[0.15, 0.075], [0.20, 0.05]);
ax2 = zeros(1,nM);
for i=1:nM

    pval_mat = pval_mat_all(:,:,i); 
    mw = mw_all(:,:,i);
    sw = sw_all(:,:,i);
    
    % Prepare P value matrix
    pval_mat(pval_mat==0)=realmin; M = -log10(pval_mat);
    Mx = M;
    Mx ( M < T(1) ) = NaN;
    for j = 1:nT-1
        Mx ( M>= T(j) & M<T(j+1) ) = j;
    end
    Mxx = nT +2.*eye(size(Mx,1)); Mxx(Mxx==nT)=NaN;
    Ixx = itril(size(Mx,1),-1);
    Jxx = itriu(size(Mx,1),1);
    Mxx(Jxx) = Mx(Jxx);  
    Mxx(Ixx) = nT+1;
    Mxx(Mxx > nT+2) = nT-1;
    mi = nM+i;
    li = i;
    
    % Print pairwise P value matrix
    imagesc(ax(mi),Mxx, 'AlphaData', ~isnan(Mxx));
    
    %ax.Position = [0.15 0.15 0.775 0.55];
    cmap1 = jet(nT-1);
    cmap2 = flipud(gray(3));
    cmap = [cmap1;cmap2];
    colormap(ax(mi),cmap);
    ax(mi).CLim =[ 1 nT+3 ];
    
    if i==nM
        CLabels = cellstr(num2str(T','%1.2f')); CLabels{end} = sprintf('>%1.2f',T(end-1));
        cbar1 = colorbar(ax(mi),'Ticks', 1:nT, 'TickLabels', CLabels); cbar1.Label.String = '-log10(P_{FDR} value)'; cbar1.Limits=[1 nT];
    end
    
    if exist('ticklabels','var') && ~isempty(ticklabels) && numel(ticklabels)==size(pvals,1)
        ax(mi).XTick = 1: numel(ticklabels);
        ax(mi).XTickLabel = ticklabels;
        ax(mi).XTickLabelRotation = 45;
        ax(mi).FontWeight='bold';
        ax(mi).FontSize=11;
        if i==1
            ax(mi).YTick = ax(mi).XTick;
            ax(mi).YTickLabel = ticklabels;
            ax(mi).TickLabelInterpreter='none';
        end
    end
    
    hold(ax(li),'on'); Gfl = false;
    if exist('G','var') && ~isempty(G), 
        Gfl = true;
        axes(ax(li)); 
        nG = size(G,3);
        if nG<2
            violin(G, 'facecolor', rgb('DimGrey'), 'edgecolor', 'none');
            xlims = size(G,2);
        else
            mG = size(G,2);
            Goff = (1:mG)/5; Goff = Goff - mean(Goff);
            for n=1:mG
                Xvec = (1:nG) - Goff(n)*-1;
                Gn = reshape(G(:,n,:),size(G,1),[]);
                notBoxPlot(Gn, Xvec,'jitter',0.15);    
            end
            xlims = size(G,3);
        end
        
    else
        plot(ax(li), mw, 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'MarkerEdgeColor', 'w');
        if size(sw,1) == 2
            er = errorbar(ax(li), 1:numel(mw), mw, sw(1,:), sw(2,:), 'LineStyle', 'none', 'Color', 'k'); 
            lgstr = {'Median','25%/75% percentiles'};
        else
            er = errorbar(ax(li), 1:numel(mw), mw, sw, 'LineStyle', 'none', 'Color', 'k'); 
            lgstr = {'Median','Stdev'};
        end
        xlims = size(mw,2);
    end
    xlim(ax(li),[0.5 xlims+0.5]);
    ax(li).XTickLabel = [];
    ax(li).Position(3) = ax(mi).Position(3);
    ax(li).Box = 'on';
    ax(li).YGrid = 'on';
    ax(li).FontWeight = 'bold';
    ax(li).FontSize = 11;
    if i==1, 
        ylabel(ax(li),perftitle);
        ax(li).YTickLabelMode='auto';
        if ~Gfl, legend(ax(li),lgstr); end
    end
    
    if exist('hlinepos','var') && ~isempty(hlinepos), 
        yline(ax(li),hlinepos,'--');
        if ~isempty((ax(li).Legend))
            ax(li).Legend.String(end)=[];
        
        end
    end
    ax(li).YAxis.LimitsMode='auto';
    if ~isempty((ax(li).Legend)); ax(li).Legend.Location='best'; end
 
    %cbar2 = colorbar(ax(li),'Ticks', 1:numel(T), 'TickLabels', CLabels); cbar2.Label.String = '-log10(P_{FDR} value)'; cbar2.Limits=[1 numel(T)]; cbar2.Visible='off';
end
