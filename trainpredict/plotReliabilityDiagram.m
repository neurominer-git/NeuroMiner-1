% Plots a reliability diagram given probabilty estimates in p, with
% corresponding true labels y
function plotReliabilityDiagram(p, y, legends, gridSize)

    if nargin < 4, gridSize = 2; end;

    mycolor = 'brm';
    marker = 'ox*';
    
    % We'll operate one column by one    
    for k = 1 : size(p, 2)    
        pk = p(:,k);

        [pDiscK, empProbK] = binProbabilities(pk, y, gridSize);

        plotStr = sprintf('%s%s',mycolor(k),marker(k));
        plot(pDiscK, empProbK, plotStr);
        hold on;        
    end

    plot([0:0.1:1], [0:0.1:1], 'k--', 'LineWidth', 4);
    hold off;

    legends{end+1} = 'Calibrated';
    xlabel('Classifier scores'); ylabel('Empirical probability of positive'); legend(legends);
    set(gca, 'FontSize', 16); set(get(gca,'XLabel'), 'FontSize', 16); set(get(gca,'YLabel'), 'FontSize', 16);    
