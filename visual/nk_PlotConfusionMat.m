function [h_mat, hStrings] = nk_PlotConfusionMat(mat, groupnames, ha, contigtable)

if ~exist('ha','var') || isempty(ha) || ~ishandle(ha)
    ha=gca;
end
ngroups                 = size(mat,1);
ngroupvec               = [3  5  7  9  11];
FontSizes               = [12 11 10 9  8];
ngroupsl                = sum(ngroupvec <= ngroups); if ~ngroupsl, ngroupsl = 1; end

FontSize                = FontSizes(ngroupsl);
h_mat                   = imagesc(mat);             % Create a colored plot of the matrix values
colormap(flipud(gray));                             % Change the colormap to gray (so higher values are
                                                    %   black and lower values are white)

textStrings             = num2str(mat(:),'%1.1f');     % Create strings from the matrix values
textStrings             = strtrim(cellstr(textStrings));  % Remove any space padding
[x,y]                   = meshgrid(1:ngroups);            % Create x and y coordinates for the strings
hStrings                = text(x(:),y(:),textStrings(:),...      % Plot the strings
                            'HorizontalAlignment','center');
set(ha,'CLim',[0 100]);
midValue                = mean(get(ha,'CLim'));     % Get the middle value of the color range
textColors              = repmat(mat(:) > midValue,1,3);  % Choose white or black for the
                                                    % text color of the strings so
                                                    % they can be easily seen over
                                                    % the background color
set(hStrings,{'Color'},num2cell(textColors,2), ...
    'FontSize',12, 'FontWeight', 'bold');           % Change the text colors

set(gca,'XTick',1:ngroups,...                       % Change the axes tick marks
        'XTickLabel',groupnames,...                 % and tick labels
        'YTick',1:ngroups,...
        'YTickLabel',groupnames,...
        'TickLength',[0 0],'FontWeight','bold', 'FontSize',FontSize);
    
if ~isempty(which('ytickangle')), ytickangle(90); end

xlabel('Predicted'); ha.XLabel.FontSize = 12;
ylabel('Observed'); ha.YLabel.FontSize = 12;
try
    ha_height = contigtable.Position(2)- ha.Position(2);
    ha.Position(4) = ha_height - 0.06;
    ha.Position(2) = 0.08;
catch
end