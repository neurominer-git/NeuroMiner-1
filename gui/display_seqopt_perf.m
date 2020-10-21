function R = display_seqopt_perf(NM, analind, Cind, NodeN)

R=[];
if ~exist('analind','var') || isempty(analind), errordlg('Please specify a sequential prediction analysis'); return; end
if ~strcmp(NM.analysis{analind}.params.TrainParam.SVM.prog,'SEQOPT'), errordlg('The specified analysis is not a sequential predictor!'); return; end
if size(NM.analysis{analind}.params.TrainParam.SVM.SEQOPT.C,1)>1 && (~exist('Cind','var') || isempty(Cind)), errordlg('The sequential predictor contains multiple sequences. Choose one!'); return; end

sz = get(0,'ScreenSize');
win_wdth = sz(3)/1.1; win_hght = sz(4)/1.2; win_x = sz(3)/2 - win_wdth/2; win_y = sz(4)/2 - win_hght/2;
ax_wdth = 0.265; ax_hght = 0.425;

% Group names
G1n = NM.groupnames{1}; G2n = NM.groupnames{2};

% Labels
L = NM.label; L(L==2) = -1;

% Prediction trajectories
P = NM.analysis{analind}.GDdims{1}.CV2grid.decvaltraj;

% Propagation thresholds
mU =  NM.analysis{analind}.GDdims{1}.grid.mean_SeqPercUpper;
sU =  NM.analysis{analind}.GDdims{1}.grid.se_SeqPercUpper;
mL =  NM.analysis{analind}.GDdims{1}.grid.mean_SeqPercLower;
sL =  NM.analysis{analind}.GDdims{1}.grid.se_SeqPercLower;

% Case propagations
I = NM.analysis{analind}.GDdims{1}.CV2grid.caseprop_node;

% Node names of sequence
if ~exist('NodeN','var') || isempty(NodeN)
    FeatN = NM.analysis{analind}.params.TrainParam.STACKING.featnames;
    NodeN = FeatN(Cind,:);
end

% Graphical elems
% Colors
cl1 = rgb('Red'); cl2 = rgb('Blue'); cl3 = rgb('DarkGreen'); cl4 = rgb('Peru'); cl5 = rgb('SaddleBrown'); cl6 = rgb('DarkGrey');

% Line Weights
lw1 = 2; lw2 = 2;

% Marker Sizes
mk1 = 12; mk2 = 12;

nI      = numel(unique(I));
R.FPRs  = zeros(1,nI);
R.BACs  = zeros(1,nI);
R.SENSs = zeros(1,nI);
R.SPECs = zeros(1,nI);
R.pLR   = zeros(1,nI);
R.PercRemain = zeros(1,nI-1);
ii      = true(size(P,1),1);

for i=1:nI
    ix = I < i;
    ii(ix) = false;
    R.Nremain(i) = sum(ii);
    R.ObsG1(i) = sum(L(ii)==1);
    R.ObsG2(i) = sum(L(ii)==-1);
    R.PredG1(i) = sum(sign(P(ii,i))==1);
    R.PredG2(i) = sum(sign(P(ii,i))==-1);
    R.ObsG1G2rat(i) = R.ObsG1(i)*100 / sum(R.ObsG1(i)+ R.ObsG2(i));
    R.PredG1G2rat(i) = R.PredG1(i)*100 / sum(R.PredG1(i)+ R.PredG2(i));
    R.PredG1tot(i) = R.PredG1(i)*100 / size(I,1);
    R.PredG2tot(i) = R.PredG2(i)*100 / size(I,1);
    R.FPRs(i) = FPR(L,P(:,i));
    R.PPVs(i) = PPV(L,P(:,i));
    R.BACs(i) = BAC(L,P(:,i));
    R.SENSs(i) = SENSITIVITY(L,P(:,i));
    R.SPECs(i) = SPECIFICITY(L,P(:,i));
    R.pLR(i) = R.SENSs(i) / (100 - R.SPECs(i));
    R.PercPopRem(i) = (R.Nremain(1)-sum(ix))*100/R.Nremain(1);

end

for i=1:nI-1
    R.PercRemain(i) = R.Nremain(i+1)*100/R.Nremain(i); R.mean_PercRemain=mean(R.PercRemain); R.std_PercRemain=std(R.PercRemain);
end

figure('Position', [win_x win_y win_wdth win_hght]);
% Graph 1: performance measures 1
R.ax(1) = axes('Position',[ 0.05 0.075 ax_wdth ax_hght ]); hold on
R.a = plot(R.ObsG1,'MarkerFaceColor',cl1,'MarkerSize',mk1,'MarkerEdgeColor','w','Marker','o','Color',cl1,'LineWidth', lw1);
R.b = plot(R.ObsG2,'MarkerFaceColor',cl2,'MarkerSize',mk1,'MarkerEdgeColor','w','Marker','o','Color',cl2,'LineWidth' ,lw1);
R.ax(1).YAxis(1).Label.String = '# patients at sequence node';
R.ax(1).YAxis(1).Label.FontWeight = 'bold';
R.ax(1).XTick = 1 : nI;
R.ax(1).XLim = [ 0.9 nI+0.1];
R.ax(1).XTickLabel = NodeN;
R.ax(1).Box = 'on';
R.ax(1).YGrid = 'on'; 
lg = {[G1n ' cases'],[G2n ' cases']}; legend(R.ax(1), lg, 'Location','best');

% Graph 2: Sample sizes
R.ax(2) = axes('Position',[ 0.05 0.525 ax_wdth ax_hght ]); hold on
R.c = plot(R.SENSs,'MarkerFaceColor',cl1,'MarkerSize',mk1,'MarkerEdgeColor','w','Marker','o','Color',cl1,'LineWidth', lw1);
R.d = plot(R.SPECs,'MarkerFaceColor',cl2,'MarkerSize',mk1,'MarkerEdgeColor','w','Marker','o','Color',cl2,'LineWidth' ,lw1);
R.ax(2).XTick = 1 : nI;
R.ax(2).XLim = [ 0.9 nI+0.1 ];
R.ax(2).XTickLabel = [];
R.ax(2).YAxis(1).Label.String = 'Sensitivity & Specificity [%]';
R.ax(2).YAxis(1).Label.FontWeight = 'bold';
R.ax(2).Box = 'on';
R.ax(2).YGrid = 'on'; 
lg = {'Sensitivity','Specificity'}; legend(R.ax(2), lg, 'Location','best');

R.ax(3) = axes('Position',[ 0.375 0.075 ax_wdth ax_hght ]); hold on
R.e = plot(R.ObsG1G2rat,'MarkerFaceColor',cl1,'MarkerSize',mk2,'MarkerEdgeColor','w','Marker','o','Color',cl1, 'LineWidth',lw2);
R.f = plot(R.PercPopRem,'MarkerFaceColor','k','MarkerSize',mk2,'MarkerEdgeColor','w','Marker','o','Color','k', 'LineWidth',lw2);
R.ax(3).YAxis(1).Label.FontWeight = 'bold';
R.ax(3).XTick = 1 : nI;
R.ax(3).XLim = [ 0.9 nI+0.1];
R.ax(3).XTickLabel = NodeN;
R.ax(3).Box = 'on';
R.ax(3).YGrid = 'on'; 
ylabel('Cases at prognostication node [%]');
R.ax(3).YAxis(1).Label.FontWeight = 'bold';
lg = {[ G1n ' patients in cases to be prognosticated [%]'],'cases to be prognosticated [%]'}; legend(R.ax(3), lg, 'Location','best');

R.ax(4) = axes('Position',[ 0.375 0.525 ax_wdth ax_hght ]); hold on
R.g = plot(R.FPRs,'MarkerFaceColor',cl3,'MarkerSize',mk2,'MarkerEdgeColor','w','Marker','o','Color',cl3, 'LineWidth',lw2);
R.h = plot(R.PPVs,'MarkerFaceColor',cl4,'MarkerSize',mk2,'MarkerEdgeColor','w','Marker','o','Color',cl4, 'LineWidth',lw2);
R.ax(4).YAxis(1).Label.FontWeight = 'bold';
R.ax(4).XTick = 1 : nI;
R.ax(4).XLim = [ 0.9 nI+0.1 ];
R.ax(4).XTickLabel = [];
R.ax(4).Box = 'on'; 
R.ax(4).YGrid = 'on'; 
ylabel('FPR & PPV [%]');
R.ax(4).YAxis(1).Label.FontWeight = 'bold';
yyaxis right;
R.f2 = plot(R.pLR,'MarkerFaceColor',cl5,'MarkerSize',mk2,'MarkerEdgeColor','w','Marker','o','Color',cl5, 'LineWidth',lw2);
R.ax(4).YAxis(2).Color = cl5;
%R.ax(4).YAxis(2).Label.String = 'Positive Likelihood Ratio';
lg = {'False Positive Rate [%]', 'Positive Predictive Value [%]','Positive Likelihood Ratio'}; legend(R.ax(4), lg, 'Location','best');

R.ax(5) = axes('Position',[ 0.7 0.075 ax_wdth ax_hght ]); hold on
R.k = plot(R.PredG1tot,'MarkerFaceColor',cl1,'MarkerSize',mk2,'MarkerEdgeColor','w','Marker','o','Color',cl1, 'LineWidth',lw2);
R.l = plot(R.PredG2tot,'MarkerFaceColor',cl2,'MarkerSize',mk2,'MarkerEdgeColor','w','Marker','o','Color',cl2, 'LineWidth',lw2);
R.ax(5).YLim = [ 0 100 ];
R.ax(5).XTick = 1 : nI;
R.ax(5).XLim = [ 0.9 nI+0.1 ];
R.ax(5).XTickLabel = NodeN;
R.ax(5).Box = 'on'; 
R.ax(5).YGrid = 'on';
ylabel('Predictions at prognostication node [%]');
R.ax(5).YAxis(1).Label.FontWeight = 'bold';
lg = {[G1n ' predictions'],[G2n ' predictions']}; legend(R.ax(5), lg, 'Location','best');

R.ax(6) = axes('Position',[ 0.7 0.55 ax_wdth ax_hght-0.025 ]); hold on
%R.i(1) = errorbar(mU, sU, 'LineWidth', lw1);
R.i(1) = plot(mU,'MarkerFaceColor','k','MarkerSize',mk2,'MarkerEdgeColor','w','Marker','o','Color','k', 'LineWidth',lw2);
%R.j(1) = errorbar(mL, sL, 'LineWidth', lw1);
R.j(1) = plot(mL,'MarkerFaceColor','k','MarkerSize',mk2,'MarkerEdgeColor','w','Marker','o','Color','k', 'LineWidth',lw2);
R.ax(6).YLim = [ 0 100 ];
R.ax(6).XTick = 1 : nI-1;
R.ax(6).XLim = [ 0.9 nI-1+0.1 ];
R.ax(6).Box = 'on'; 
R.ax(6).YGrid = 'on';
ylabel('Propagation thresholds [percentiles]');
R.ax(6).YAxis(1).Label.FontWeight = 'bold';
for i=1:nI-1, NodeNx{i} = sprintf('Node %g=>%g', i, i+1); end
R.ax(6).XTickLabel = NodeNx;
yline(50,'k:','LineWidth',lw2)

for i=1:numel(R.ax)
   R.ax(i).Color = [0.95 0.95 0.95]; 
end
