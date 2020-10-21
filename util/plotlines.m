function plotlines(Data, Label, LabelNames, Featnames)

figure; hold on
[n,m,o]=size(Data);
uL = unique(Label);
nuL = numel(uL);
colpt = {'b-','r-','g-','m-','c-','k-'};
colpto = {'bo-','ro-','go-','mo-','co-','ko-'};
for j=1:m
    if o>1, 
        ax_j=subplot(1,m,j); hold on; 
    else
        ax_j = axes;
    end
    if exist('Featnames','var') && ~isempty(Featnames), title(Featnames{j}); end
    for i=1:n
        if o>1
            iY = Data(i,j,:);
        else
            iY = Data(i,:);
        end
        colpt_i = colpt{Label(i)};
        p=plot(ax_j, 1:numel(iY(:)), iY(:), colpt_i,'LineWidth',0.05,'MarkerSize',5);
        p.Color(4)=0.3;
    end

    for i=1:nuL
        if o>1
            mY = nanmean(Data(Label==uL(i),j,:));
            sY = nanstd(Data(Label==uL(i),j,:));
        else
            mY = nanmean(Data(Label==uL(i),:));
            sY = nanstd(Data(Label==uL(i),j,:));
        end
        colpt_i = colpto{Label(i)};
        p(i)= plot(ax_j, 1:numel(mY(:)), mY(:), colpt_i,'LineWidth',2,'MarkerSize',12,'MarkerFaceColor', colpt_i(1) );
    end
    if j==1; legend(p,LabelNames); end
end

