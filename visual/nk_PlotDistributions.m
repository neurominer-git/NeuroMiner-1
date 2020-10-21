function [fg, h, hp] = nk_PlotDistributions(NM, analind, singleP)

P = NM.analysis{analind}.GDdims{1}.BinClass{1}.mean_predictions;
L = NM.label;
vec=1:numel(unique(L));
fg=figure;hold on
h=cell(1,numel(vec)); hp = h;
cl = {'r','b','g','c','m'};
xline(0,'LineWidth',1.5,'Color',[0.6 0.6 0.6],'LineStyle',':')
for i=1:numel(vec)
    h{i} = histfit(P(L==vec(i)),10,'kernel');
    delete(h{i}(1));
    if numel(vec)==2 && i==2
        h{i}(2).YData = -1*h{i}(2).YData*100/sum(L==vec(i));
    else
        h{i}(2).YData = h{i}(2).YData*100/sum(L==vec(i));
    end
    h{i}(2).Color = cl{i};
    h{i}(2).Color(4) = 0.3;
    hp{i} = patch([min(h{i}(2).XData)-0.00001,h{i}(2).XData,max(h{i}(2).XData)+0.00001],[0,h{i}(2).YData,0],cl{i},'FaceAlpha',0.1,'EdgeColor','none');
end
ax= gca;
ax.XAxis.Label.String='Decision score';
ax.YAxis.Label.String='Percentage of outcome class';
if exist('singleP','var') && ~isempty(singleP)
    er= errorbar(ax,singleP(1),0,singleP(2),singleP(3),'ko-');
    er.MarkerSize=12;
    er.MarkerFaceColor='k';
end
view([90 -90])