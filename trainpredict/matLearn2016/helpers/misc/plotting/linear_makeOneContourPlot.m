function linear_makeOneContourPlot(Xtrain,ytrain,model)
[n,d] = size(Xtrain);

% adjust dimension for bias
d = d-1;
if d == 2
    increment = 100;
    
    nClasses = model.nClasses;
    name = model.name;
    
    % Plot Data
    figure;
    clf; hold on;
    colors = getColorsRGB;
    for c = 1:nClasses
        plot(Xtrain(ytrain==c,2),Xtrain(ytrain==c,3),'.','color',...
            colors(c,:));
    end
    
    domain1 = xlim;
    domain1 = domain1(1):(domain1(2)-domain1(1))/increment:domain1(2);
    domain2 = ylim;
    domain2 = domain2(1):(domain2(2)-domain2(1))/increment:domain2(2);
    
    d1 = repmat(domain1',[1 length(domain1)]);
    d2 = repmat(domain2,[length(domain2) 1]);
    
    
    %fprintf('Plotting %s...\n',model.name);
    
    % Compute values
    vals = model.predict(model,[ones(numel(d1),1) d1(:) d2(:)]);
    % figure;clf;hold on;
    
    zData = reshape(vals,size(d1));
    u = unique(zData(:));
    
    % Remove classes that don't occur (for plotting purposes)
    for c = nClasses:-1:1
        if ~any(zData(:)==c)
            zData(zData > c) = zData(zData > c) - 1;
        end
    end
    
    % Set Colormap
    cm = colors(u,:)/2;
    colormap(cm);
    % Make Contour Plot
    contourf(d1,d2,zData,1:max(zData(:)),'k');
    for c = 1:nClasses
        plot(Xtrain(ytrain==c,2),Xtrain(ytrain==c,3),'.','color',colors(c,:));
    end
    
    hold off;
    % Make Legend  
    %legend(name{m});
    title(sprintf('%s (train err = %.3f)',model.name,model.trainError));
    %fprintf('Press any key to continue')
end
end