function [] = plotRegression1D(X,y,varargin)
% plotRegression1D(X,y,model1,model2,...)
%
% Displays the 1D regression problem data as well the predictions of different
% models
[nTrain,nFeatures] = size(X);

assert(nFeatures==1,'plotRegression1D can only plot 1D data sets');

nIncrements = 100;


figure;
h(1,1) = plot(X,y,'k.');
names{1,1} = 'Train';
hold on;

xRange = xlim;
xRange = [xRange(1):(xRange(2)-xRange(1))/nIncrements:xRange(2)]';

colors = getColorsRGB;
for m = 1:length(varargin)
    model = varargin{m};
    yhat = model.predict(model,xRange);
    h(m+1,1) = plot(xRange,yhat,'LineWidth',3,'Color',colors(m,:));
    names{m+1,1} = model.name;
end
legend(h,names);

