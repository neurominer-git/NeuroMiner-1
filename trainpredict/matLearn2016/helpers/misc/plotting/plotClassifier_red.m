function plotClassifier_2(X,y,model)

[n,p] = size(X);
nClasses = max(y) + 1;
if p ~= 2
    fprintf('Plotting only supported for 2D data\n');
    return;
end

hold on
colors = getColorsRGB;
for c = 1:nClasses
    plot(X(y==c,1),X(y==c,2),'.','color',colors(c,:));
end

increment = 50;
domainx = xlim;
domain1 = domainx(1):(domainx(2)-domainx(1))/increment:domainx(2);
domainy = ylim;
domain2 = domainy(1):(domainy(2)-domainy(1))/increment:domainy(2);
d1 = repmat(domain1',[1 length(domain1)]);
d2 = repmat(domain2,[length(domain2) 1]);

yhat = model.predict(model, [d1(:) d2(:)]);

z = reshape(yhat,size(d1));
u = unique(z(:));
% For plotting purposes, remove classes that don't occur
for c = nClasses:-1:1
    if ~any(z(:)==c)
        z(z > c) = z(z > c)-1;
    end
end
cm = colors(u,:)/2;
colormap(cm);
contourf(d1,d2,z,1:max(z(:)),'r');

colors = getColorsRGB;
for c = 1:nClasses
    plot(X(y==c,1),X(y==c,2),'.','color',colors(c,:));
end

xlim(domainx);
ylim(domainy);
%legend({'Class 1','Class 2',method});
title(model.name);