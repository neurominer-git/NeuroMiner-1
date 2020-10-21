% plots for BIBM 2012

% sparse coding
dataSets={'Breast';'Prostate'};
classifiers={'KNNLS';'Kl_1NNLS';'Kl_1LS';'NNLS';'l_1NNLS';'l_1LS';'l_1LS-IP';'l_1LS-PX';'kNN';'SVM'};
Acc=zeros(2,10);
%          (,2^-1,)(2^-1,2^-10)(2^0,2^-9)       2^-9    (,,2^-10)                  k=1                                              
Acc(1,:)= [0.7865,  0.7862,  0.7915,  0.7804,  0.7804,  0.7822,  0.7833,  0.7896,  0.7553,  0.7798]; % linearSVM 0.8072

%          2^3     (2^0,2^-6)                   2^-10    2^-10                      k=3                            
Acc(2,:)= [0.9272,  0.9264,  0.9545,  0.9319,  0.9317,  0.9666,  0.9689,  0.8568,  0.8780,  0.9274]; % linearSVM 0.9398

std=zeros(2,10);
std(1,:)= [0.0102,  0.0147,  0.0137,  0.0128,  0.0131,  0.0173,  0.0168,  0.0109,  0.0251,  0.0187]; % 0.0224
std(2,:)= [0.0099,  0.0100,  0.0077,  0.0155,  0.0070,  0.0100,  0.0057,  0.0053,  0.0108,  0.0087]; % 0.0089

time=zeros(2,10);
time(1,:)=[0.4598,  0.7591,  1.1373,  0.5191,  0.5206,  0.9025,  38.4890, 0.6517,  0.1051,  3.6553]; % 3.1297 % check the first three time data
time(2,:)=[3.7570,  3.8729,  37.8136, 3.7346,  3.8445,  52.5262, 705.4383,2.6389,  0.4286,  8.6356];% 5.4493

% reorder
reOrder=[6,7,8,4,5,3,1,2,9,10];
classifiers=classifiers(reOrder);
Acc=Acc(:,reOrder);
std=std(:,reOrder);
time=time(:,reOrder);

xDescription='Data';
yDescription='Accuracy';
saveFigName='BreastProstateSC.eps';
plotBarError(Acc,std,dataSets,classifiers,xDescription,yDescription,saveFigName);
saveFigName='BreastProstateSCTime.eps';
plotTime(time(2,:),classifiers,saveFigName);
xDescription='Data';
yDescription='log_2(Computing Time)';
plotBarError(log2(time),[],dataSets,classifiers,xDescription,yDescription,saveFigName);
save('BreastProstate.mat','Acc','std','time','dataSets','classifiers');
%% SR
dataSets={'Breast';'Prostate'};
classifiers={'KSR-NNLS';'KSR-l_1NNLS';'KSR-l_1LS';'SR-NNLS';'SR-l_1NNLS';'SR-l_1LS';'semi-NMF';'SR-l_1LS-IP';'SR-l_1LS-PX'};
Acc=zeros(2,9);
%          (4,2^4,)(4,2^0,-6)(4,2^0,-6) (4)   (4,,2^-8) (4,,2^-6)
Acc(1,:)= [0.8361,  0.8373,  0.8367,  0.8323,  0.8341,  0.8392,  0.8346,  0.8311,  0.2970];
%          (24,2^2)(26,2^1,2^-11)(22,2^1,2^-11)(27)(25,,2^-11)(24,,2^-8)                                       
Acc(2,:) =[0.9290,  0.9280,  0.9277,  0.9295,  0.9312,  0.9342,  0.9292,  0.8934,  0.7139]; % prostate

std=zeros(2,9);
std(1,:)= [0.0199,  0.0157,  0.0184,  0.0178,  0.0155,  0.0182,  0.0153,  0.0204,  0.0516];
std(2,:)= [0.0079,  0.0087,  0.0091,  0.0082,  0.0084,  0.0077,  0.0085,  0.0383,  0.0799];
time=zeros(2,9);
time(1,:)=[0.0918,  0.1498,  0.6605,  0.0909,  0.0914,  0.7172, 31.4393,  69.8150, 9.3556];
time(2,:)=[1.9306,  2.4136,  5.9025,  2.1299,  3.4172,  13.3478,256.1137, 68.1992, 15.9468];

% reorder
reOrder=[6,8,9,4,7,5,3,1,2];
classifiers=classifiers(reOrder);
Acc=Acc(:,reOrder);
std=std(:,reOrder);
time=time(:,reOrder);

xDescription='Data';
yDescription='Accuracy';
saveFigName='BreastProstateSR.eps';
plotBarError(Acc,std,dataSets,classifiers,xDescription,yDescription,saveFigName);
saveFigName='BreastProstateSRTime.eps';
plotTime(time(1,:),classifiers,saveFigName);
xDescription='Data';
yDescription='log_2(Computing Time)';
plotBarError(log2(time),[],dataSets,classifiers,xDescription,yDescription,saveFigName);
save('BreastProstateSR.mat','Acc','std','time','dataSets','classifiers');