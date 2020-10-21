function [] = plot2DLine(model)

increment = 500;

%figure;clf;
domain1 = xlim;
domain1 = domain1(1):(domain1(2)-domain1(1))/increment:domain1(2);
domain2 = ylim;
domain2 = domain2(1):(domain2(2)-domain2(1))/increment:domain2(2);

d1 = repmat(domain1',[1 length(domain1)]);
d2 = repmat(domain2,[length(domain2) 1]);

vals = model.predict(model,[d1(:) d2(:)]);


zData = reshape(vals,size(d1));
contourf(d1,d2,zData+rand(size(zData))/1000,[-1 0],'k');