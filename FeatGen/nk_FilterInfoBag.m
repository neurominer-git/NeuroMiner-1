function [MinInfY, symY, binY] = nk_MIBag(Y, labels, SymFlag)
%
% function [MinInfY, symY, binY] = nk_FilterInfoBag(Y, labels, SymFlag, InfoType, nperms)
%
% Inputs:
% -------
% Y:                Data (rows = observations, columns = features)
% labels:           Class memberships
% SymFlag = 0/1:    Data is symbolized first
% _________________________________
% (c) Nikolaos Koutsouleris 05/2009

if nargin < 4
    InfoType = 1;
end

if nargin < 3
	SymFlag = 1;
end

kX = size(Y,1);ix = size(Y,2);
fprintf('\n% g oberservations, %g features',kX,ix)
MinInfY = zeros(ix,1);

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

try 
    matlabpool open
    try
        parfor i=1:ix % perform parallel computation

            sY = symY(:,i);
            lb = labels;

            if any(symY(:,i))

                BagY = zeros(1,kX);

                % Perform generalization analysis using LOO
                for h=1:kX
                    Ind = 1:kX; Ind(h) = [];
                    BagY(h) = mutualinfo(sY(Ind),lb(Ind));

                end
                MinInfY(i) = min(BagY);% use the minimum value
            else
                MinInfY(i) = 0;
            end
        end
        flg=1;
    catch
        flg=1;
    end
    if flg, matlabpool close; end

catch

    for i=1:ix

        sY = symY(:,i);

        if any(symY(:,i))

            BagY = zeros(1,kX);

            % Perform generalization analysis using LOO
            for h=1:kX
                Ind = 1:kX; Ind(h) = [];
                BagY(h) = mutualinfo(sY(Ind),labels(Ind));

            end
            MinInfY(i) = min(BagY);% use the minimum value
        else
            MinInfY(i) = 0;
        end
    end
end
return