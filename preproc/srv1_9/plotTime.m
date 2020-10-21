function plotTime(time, classifiers,saveFigName)

figure('OuterPosition',[0 0 850 600],...
'MenuBar','figure',...
'ToolBar','figure',...
'PaperPositionMode','auto');
bar(1:numel(time),log2(time),'b'); 
set(gca,'XTickLabel',classifiers);
set(gca,'FontSize',10);
xlabel('Approach','FontSize',14);
ylabel('log_2(Computing Time)','FontSize',14);
print(gcf,'-depsc2','-r300',saveFigName);
end