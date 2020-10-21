function TestMe(No_RandomFeature)

%Conduct experiment on the well-known spiral data
%---------------------------------------------------------------|
%Input:                     
%      No_RandomFeature: number of randomly distributed features
%---------------------------------------------------------------
close all

%load the data
load spiral.mat
N = length(targets);             % Number of patterns
Original_dim = size(patterns,1); % Number of original features
targets = targets+1;             % class label {1,2}
patterns = [patterns; randn(No_RandomFeature,N)]; %Add some irrelevant features
dim = size(patterns,1);          % Data dimenionality

%Preprocess the data: tranform each feature into [0, 1] so that features are comparable. 
%This is a commonly used pre-processing step in nearest-neighbor algorithm.
[MIN,I] = min(patterns,[],2);
[MAX,I] = max(patterns,[],2);
for n=1:dim
    if (MAX(n)-MIN(n))==0
        patterns(n,:)=0;
    else
        patterns(n,:) = (patterns(n,:)-MIN(n))/(MAX(n)-MIN(n));
    end
end

%plot the first two features
index_1 = find(targets==1);
index_2 = find(targets==2);

figure(1)
plot(patterns(1,index_1),patterns(2,index_1),'*','MarkerSize',10);
hold on
plot(patterns(1,index_2),patterns(2,index_2),'ro','MarkerSize',10);
title('The first two features')
axis square;axis tight
boldify1
drawnow

%parameters
Para.plotfigure = 1;     % 1: plot of the result of each iteration; 0: do not plot
Para.distance = 'block'; % 'euclidean';  
Para.sigma= 2;           % kernel width; If the algorithm does not converge, use a larger kernel width.
Para.lambda = 1;         % regularization parameter
%I arbitarily set sigma= 2 and lambda = 1. The proposed algorithm is not sensitive to parameters. 
%The algorithm can used for classification. The parameters can be learning via cross-validation (see the paper).

s = cputime;
Weight = IMRelief_Sigmoid_FastImple(patterns, targets(:), Para);
CPUTime = cputime-s;
disp(['>>> The total CPU time is ' num2str(CPUTime) '.'])

figure;
semilogx(Weight/max(Weight),'-o','LineWidth',1,'MarkerFaceColor','w','MarkerSize',10)
hold on
plot([Original_dim,Original_dim],[0,1],'r--', 'LineWidth',2);
xlabel('Features')
ylabel('Feature Scores')
title('Feature Weights');
axis tight
boldify1

return
%=====================End of The Code=====================================



