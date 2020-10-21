function [infY, alpha, sig, infYPerm,symY, binY] = nk_FilterInfo(Y, labels, SymFlag, InfoType, nperms, alphacut)
%
% function [infY, alpha, sig, infYPerm, symY, binY] = nk_FilterInfo(Y, labels, SymFlag, InfoType, nperms)
%
% Inputs:
% -------
% Y:                Data (rows = observations, columns = features)
% labels:           Class memberships
% SymFlag = 0/1:    Data is symbolized first
% InfoType = 1/2:   Mutual Information/Conditional Entropy
% nperms:           # permutations to assess feature significance [500]
% alphacut:         cutoff for significance [0.05] = 95% CI
% _________________________________
% (c) Nikolaos Koutsouleris 05/2009

if nargin < 6
    alphacut = .05;
end

if nargin < 5
    nperms = 500; % Permutation array for feature selection
end

if nargin < 4
    InfoType = 1;
end

if nargin < 3
	SymFlag = 1;
end

kX = size(Y,1);ix = size(Y,2);
fprintf('\n% g oberservations, %g features',kX,ix)
infY = zeros(ix,1);
alpha = zeros(ix,1);	

if SymFlag
    fprintf('\nSpecify the encoding strategy:')
    seql = nk_input('Sequence length',0,'e',5);
    binstart = nk_input('Minimum bins',0,'e',3);
    binstop = nk_input('Maximum bins',0,'e',3);
    [symY, binY] = nk_Symbolize(Y, seql, binstart, binstop);
else
	symY = Y;
end

% Compute conditional entropy of each symbolized feature and the class
% labels for the training and test data

switch InfoType
    
    % Mutual Information is used
    case 1 
        
        for i=1:ix  
            % The first column is the original observation
            if any(symY(:,i))
                
                infY(i) = mutualinfo(symY(:,i),labels);

                PermY = zeros(1,nperms);
                for h=1:nperms % perform permutation analysis
                    p1Ind = randperm(kX);
		    p2Ind = randperm(kX);
                    PermY(h) = mutualinfo(symY(p1Ind,i),labels(p2Ind));
                end
                alpha(i) = sum(PermY>infY(i))/nperms;
                fprintf('\nComputing mutual information of feature %g: a=%1.2f',i,alpha(i))
            else
                infY(i) = 0;
                alpha(i)=1;
            end
            
        end
        %mRMRTrain = mrmr_mid_d(symY,labels,ix);
        [sinfY, sI] = sort(alpha,'ascend');
        sig = alpha <= alphacut;
        
    % Conditional entropy is used    
    case 2 
        for i=1:ix

            fprintf('\nComputing conditional entropy of feature %g',i)
            if any(symY(:,i))

                infY(i) = condentropy(labels,symY(:,i));
            else
                infY(i) = 0;
            end
        end
        [sinfY, sI] = sort(infY,'ascend');
end

figure(2)
subplot(3,1,1); plot(infY,'b.')
subplot(3,1,2); plot(sinfY,'b.')
subplot(3,1,3); plot(alpha(sI),'b.')
return