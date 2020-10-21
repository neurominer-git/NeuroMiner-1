function [classes, perf, dist, sim, minimum] = nk_MultiDecideErrorCorrOutCodes(X, L, ClX, Groups, Coding, Decoding, WeightFlag, ProbComp)
% [classes, ECOC] = nk_ErrorCorrOutCodes(X, L, ClX, Coding, Decoding, Weightflag)
% ===================================================================================
% 
% This function uses error-correcting outp codes to compute multi-class
% group memberships for data instances in X. X contains the concatenated
% output functions of an ensemble of binary classifiers
% 
% Inputs:
% =======
% X         : The data instances to be classified, where rows (n) are instances
%             and columns (m) the decision outputs of the base learners
% L         : Label vector
% ClX       : A 1 x m vector identifying the dichotomization index of the
%             respective base learners
% Coding    : Coding strategy. Currently: 
%               1 = one-vs-one 
%               2 = one-vs-all 
% Decoding  : Decoding strategy. Currently: 
%               1 = Hamming distance 
%               2 = Euclidean Decoding
%               3 = Laplacian Decoding
%               4 = Attenuated Euclidean
%               5 = Linear Loss-based Decoding
% WeightFlag: Optionally, weight distances with number of dichotomizers in
%             the respective class
% 
% Outputs:
% ========
% classes   : predicted class-membership for X
% perf      : prediction performance
% dist      : distances for each class
% minimum   : minimum distance vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 02/2011

%% MAIN ECOC routine
% Get # of classes
global SVM 

if ~isempty(SVM) && SVM.GridParam == 14
        multimode = 1;
    else
        multimode = 0;
end

if ~exist('ProbComp','var') || isempty(ProbComp)
    ProbComp = 'invnormquad';
end

nr_classes = numel(unique(ClX));

% Only binary decoding needed
if nr_classes==1, nr_classes = 2; end;

% Compute code words if not available
switch Coding
    case 1 % One-vs-One
        %fprintf('\nOne-vs-One coding')
        oECOC = nk_OneVsOne(ClX, Groups);

    case 2 % One-vs-All
        oECOC = nk_OneVsAll(ClX, Groups);
end

% if WeightFlag
%     Weights = zeros(1,nr_classes);
%     for curclass=1:nr_classes
%         Weights(curclass) = sum(ClX==curclass);
%     end
%     [~,mXI] = max(Weights);
%     Weights = Weights./Weights(mXI);
%     
% end

% Decode (with or without weighting according to number of dichotomizers in each class)
[classes, dist, minimum] = decoding_main(X, oECOC, Decoding);

if isempty(L), L = ones(numel(classes),1); end

switch ProbComp
    case 'softmax'
        sim = 1-nk_softmax(dist);
    case 'invnorm'
        sim = (1./dist)./(sum(1./dist,2));
    case 'invnormquad'
        sim = (1./dist.^2)./(sum(1./dist.^2,2));
end
% Compute multi-class accuracy
perf = nk_MultiPerfQuant(L,classes,multimode);

return

%% DECODING MAIN
function [classes, dist, minimum] = decoding_main(X, ECOC, decoding)

n_test = size(X,1);
n_classes = size(ECOC,1);
%classes = zeros(n_test,1);
dist = zeros(n_test, n_classes);

% Get code word for test instace X_i
for k=1:n_classes
    %Get code word from ECOC for current class
    y=ECOC(k,:); Y=repmat(y,n_test,1);
    
    switch decoding
        case 1 % Hamming distance
            dist(:,k) = HD(X,Y);
        case 2 % Euclidean distance
            dist(:,k) = ED(X,Y);
        case 3 % Laplacian decoding
            dist(:,k) = LAP(X,Y);
        case 4 % Attenuated euclidean decoding
            dist(:,k) = AED(X,Y);
        case 5 % Linear loss-based decoding
            dist(:,k) = LLB(X,Y);
%         case 6 % Cosine distance-based decoding
%             dist(:,k) = COSD(X,Y);
    end
end

[minimum, classes] = min(dist,[],2);
classes(isnan(minimum))=NaN;
return

%% DECODING FUNCTIONS
% Hamming distance decoding
function d=HD(x,y)
d=sum(abs((x-y)),2)/2;
return

% Euclidean distance decoding
function d=ED(x,y)
d=sqrt(sum((x-y).^2,2));
return

% Laplacian decoding
function d=LAP(x,y)
C=sum(x==y,2);
E=sum(x~=y,2)-sum(y==0,2);
K=2;
d=(C+E+K)./(C+1);
return

% Attenuated Euclidean decoding
function d=AED(x,y)
d=sqrt(sum((x-y).*abs(y).^2,2));       
return

% Linear Loss-based decoding
function d=LLB(x,y)
d=sum(-1*(x.*y),2);            
return

% % Cosine-distance based decoding
% function d=COSD(x,y)
% d=pdist2(x',y','cosine');
% return