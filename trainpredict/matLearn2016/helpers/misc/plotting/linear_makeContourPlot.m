
if p == 2
    increment = 100;

    % Plot Data
    figure(1);
    clf; hold on;
    colors = getColorsRGB;
    for c = 1:nClasses
        plot(X(y==c,2),X(y==c,3),'.','color',colors(c,:));
        if doTest
            plot(Xtest(ytest==c,2),Xtest(ytest==c,3),'x','color',colors(c,:));
        end
    end

    domain1 = xlim;
    domain1 = domain1(1):(domain1(2)-domain1(1))/increment:domain1(2);
    domain2 = ylim;
    domain2 = domain2(1):(domain2(2)-domain2(1))/increment:domain2(2);

    d1 = repmat(domain1',[1 length(domain1)]);
    d2 = repmat(domain2,[length(domain2) 1]);

    % Plot results
    for m = 1:length(name)
        fprintf('Plotting %s...\n',name{m});

        % Compute values
        vals = model{m}.predictFunc(model{m},[ones(numel(d1),1) d1(:) d2(:)]);


            figure(m);clf;hold on;

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
                plot(X(y==c,2),X(y==c,3),'.','color',colors(c,:));
                if doTest
                    plot(Xtest(ytest==c,2),Xtest(ytest==c,3),'x','color',colors(c,:));
                end
            end

            % Make Legend

            %legend(name{m});
            
            title(sprintf('%s (err = %.3f)',name{m},errtest(m)));
        pause;
    end


end