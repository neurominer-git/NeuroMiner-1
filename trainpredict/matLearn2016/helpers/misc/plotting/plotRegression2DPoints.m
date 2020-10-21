function [] = plotRegression2DPoints(X,y,varargin)

xlin = linspace(min(X(:,1)),max(X(:,1)),20);
ylin = linspace(min(X(:,2)),max(X(:,2)),20);
xGrid = [];
for i = 1:length(xlin);
    for j = 1:length(ylin);
        xGrid = [xGrid; xlin(i) ylin(j)];
    end
end

for m = 1:length(varargin)
    model = varargin{m};
    h(1,1) = plot3(X(:,1),X(:,2),y,'.r','MarkerSize',15); 
    yhat = model.predict(model,xGrid);
        names{1,1} = 'Train';
    hold on
    h(m+1,1) = scatter3(xGrid(:,1),xGrid(:,2),yhat,40,yhat,'filled');
    axis tight; hold off;

    names{m+1,1} = model.name;
    hold off
       legend([h(1,1), h(m+1,1)], names{1,1},'model');
title(names{m+1,1});
end

end
