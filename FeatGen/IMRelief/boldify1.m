function boldify1(LineWidth)
%BOLDIFY        Make lines and text bold.
%               BOLDIFY boldifies the lines and text of the current figure.
%               BOLDIFY(H) applies to the graphics handle H.

%               S. T. Smith

if nargin < 1, LineWidth = 2; end
h = gcf;

ha = get(h,'Children');

% Get the paper position of the current figure
units = get(gcf,'PaperUnits');
set(gcf,'PaperUnits','inches');
p_inch = get(gcf,'PaperPosition');
set(gcf,'PaperUnits',units);

for i=1:length(ha)
  hc = get(ha(i),'Children');
  for j=1:length(hc)
    chtype = get(hc(j),'Type');
   if chtype(1:4) == 'text'
     set(hc(j),'FontSize',14);         % 14 pt descriptive labels
%     set(hc(j),'FontWeight','Bold');
    elseif chtype(1:4) == 'line'
      set(hc(j),'LineWidth',LineWidth); 
    end
  end
  if get(ha(i),'Type') == 'axes'        % Axis format
    % Determine the scale in inches
    units = get(ha(i),'Units');         % Units setting on axis
    set(ha(i),'Units','normalized');
    p_norm = get(ha(i),'Position');     % Normalized Position
    set(ha(i),'Units',units);           % Back to original setting
    scale = 1/(p_norm(3)*p_inch(3));    % scale*inches -> normalized units

    set(ha(i),'FontSize',14);           % Tick mark and frame format
%    set(ha(i),'FontWeight','Bold');
    set(ha(i),'LineWidth',1);
    set(ha(i),'TickLength',[1/8 2.5*1/8]*scale); % Gives 1/8" ticks

    set(get(ha(i),'XLabel'),'FontSize',16);
%    set(get(ha(i),'XLabel'),'FontWeight','Bold');
    set(get(ha(i),'XLabel'),'VerticalAlignment','top');

    set(get(ha(i),'YLabel'),'FontSize',16);
%    set(get(ha(i),'YLabel'),'FontWeight','Bold');
    set(get(ha(i),'YLabel'),'VerticalAlignment','baseline');

    set(get(ha(i),'Title'),'FontSize',16);
%    set(get(ha(i),'Title'),'FontWeight','Bold');
  end
end