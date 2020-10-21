function D = nk_PrintCV2BarChart2(D, multiflag)

WinTag = 'PrintCVBarsBin';

if isfield(D,'h'), hx = ishandle(D.h); else hx = 0; end
if ~hx  
    h = findobj('Tag',WinTag);
    if isempty(h)
        D.h = figure('Name', D.binwintitle, ...
            'NumberTitle','off', ...
            'Tag',WinTag, ...
            'MenuBar','none', ...
            'Position', D.figuresz, ...
            'Color', [0.9 0.9 0.9]);
    else
        D.h = h; set(0,'CurrentFigure',h); clf;
    end
else
    h = D.h; set(0,'CurrentFigure',h)
end

fig = gcf; 
fig.Name = D.binwintitle;

if ~isfield(D,'hl') || ~ishandle(D.hl)
    hl = findobj('Tag','CurrParam');
    if isempty(hl)
        D.hl = axes('Parent',D.h,'Position',[0.1 0.95 0.85 0.025],'Tag','CurrParam', 'Visible','on',  'YTick', []); 
    else
        D.hl = hl;
    end
end
hold on
set(h,'CurrentAxes',D.hl); cla; barh(D.hl,1,D.pltperc, 'FaceColor', 'b'); hold on

if ~isempty(D.Pdesc{1}),
    tx = sprintf('%s\nParams: ',D.s); 
    txi = []; i=1;
    for j = 1:numel(D.Pdesc{i})
        if iscell(D.P{i}(j)), ParVal = D.P{i}{j}; else, ParVal = D.P{i}(j); end
        if isnumeric(ParVal), ParVal = num2str(ParVal,'%g'); end
        txi = [ txi sprintf('%s = %s, ', D.Pdesc{i}{j}, ParVal) ];
    end
    tx = [tx txi(1:end-2)];
else
    tx = D.s;
end

xlabel(D.hl, tx); xlim([0 100]); hold off

for p=1:numel(D.ax)
    
    % Create axes if necessary
    if ~isfield(D.ax{p},'h') || ~ishandle(D.ax{p}.h)
        ha = findobj('Tag',D.ax{p}.title);
        if isempty(ha)
            D.ax{p}.h = axes('Parent',D.h, ...
                'Position',D.ax{p}.position);
        else
            D.ax{p}.h = ha; set(h,'CurrentAxes',D.ax{p}.h);
        end
        cla;    
        % Plot bar & error charts
        [D.ax{p}.h_bar, D.ax{p}.err_bar] = barwitherr(D.ax{p}.std_y, D.ax{p}.val_y); 
        hold on
        set(D.ax{p}.h, 'Tag',D.ax{p}.title, 'FontWeight','bold'); box on
        if D.nclass ==1, set(D.ax{p}.h_bar,'FaceColor',D.ax{p}.fc); end

        % Set Y and X axis labels
        ylabel(D.ax{p}.h,D.ax{p}.ylb);

        % Set Y and X axis scaling
        if isempty(D.ax{p}.ylm)
            ylm(1) = min(D.ax{p}.val_y(:)); ylm(2) = max(D.ax{p}.val_y(:));
        else
            ylm = D.ax{p}.ylm;
        end

        ylim(ylm);
        if D.nclass > 1 && size(D.ax{p}.val_y,1) == 2
            xmax = 1;
        else
            xmax = size(D.ax{p}.val_y,1);
        end
        xlim([0.5 xmax+0.5]);
        set(D.ax{p}.h,'YTickMode','auto');
        set(D.ax{p}.h,'XTick',1:1:numel(D.ax{p}.label));
        set(D.ax{p}.h,'XTickLabel',D.ax{p}.label);
        title(D.ax{p}.h, D.ax{p}.title)
        if ~isempty(D.ax{p}.lg), 
            h_legend = legend(D.ax{p}.h, D.ax{p}.lg); 
            legend(D.ax{p}.h, 'boxoff');
            set(h_legend,'FontSize',8);
        end
        hold off
    else
        set(h,'CurrentAxes',D.ax{p}.h);
        for i = 1:D.nclass
            if sum(any(isnan(D.ax{p}.val_y(:,i)))),continue; end
            if D.nclass==1, 
                std_y = D.ax{p}.std_y';
            else
                std_y = D.ax{p}.std_y(:,i);
            end
            set(D.ax{p}.h_bar(i),'YData',D.ax{p}.val_y(:,i));
            set(D.ax{p}.err_bar(i),'YData',D.ax{p}.val_y(:,i));
            if std_y 
                set(D.ax{p}.err_bar(i),'LData',std_y);
                set(D.ax{p}.err_bar(i),'UData',std_y);
            end
        end
    end
    
    
end

if multiflag
    
    for p=1:numel(D.m_ax)
    
        % Create axes if necessary
        if ~isfield(D.m_ax{p},'h') || ~ishandle(D.m_ax{p}.h)
            hm = findobj('Tag',D.m_ax{p}.title);
            if isempty(hm)
                D.m_ax{p}.h = axes('Parent',D.h, ...
                    'Position',D.m_ax{p}.position, ...
                    'Tag',D.m_ax{p}.title, ...
                    'FontWeight','bold') ; box on
                %title(D.m_ax{p}.h, D.m_ax{p}.title)
            else
                D.m_ax{p}.h = hm;
            end
            %axes(D.m_ax{p}.h); 
            cla; hold on
            % Plot bar & error charts
            [D.m_ax{p}.h_bar, D.m_ax{p}.err_bar] = barwitherr(D.m_ax{p}.std_y, D.m_ax{p}.val_y); 
            if D.nclass ==1, set(D.m_ax{p}.h_bar,'FaceColor',D.m_ax{p}.fc); end

            % Set Y and X axis labels
            ylabel(D.m_ax{p}.h,D.m_ax{p}.ylb);

            % Set Y and X axis scaling
            if isempty(D.m_ax{p}.ylm)
                ylm(1) = min(D.m_ax{p}.val_y(:)); ylm(2) = max(D.m_ax{p}.val_y(:));
            else
                ylm = D.m_ax{p}.ylm;
            end

            ylim(ylm);
            xlim([0.5 length(D.m_ax{p}.val_y)+0.5]);
            set(D.m_ax{p}.h,'YTickMode','auto');
            set(D.m_ax{p}.h,'XTick',1:1:numel(D.m_ax{p}.label));
            set(D.m_ax{p}.h,'XTickLabel',D.m_ax{p}.label);
            title(D.m_ax{p}.h, D.m_ax{p}.title)
            hold off
        else
            set(h,'CurrentAxes',D.m_ax{p}.h);
            set(D.m_ax{p}.h_bar,'YData',D.m_ax{p}.val_y); 
            set(D.m_ax{p}.err_bar,'YData',D.m_ax{p}.val_y);
            set(D.m_ax{p}.err_bar,'LData',D.m_ax{p}.std_y');
            set(D.m_ax{p}.err_bar,'UData',D.m_ax{p}.std_y');
        end
        
    end
end 
drawnow
end