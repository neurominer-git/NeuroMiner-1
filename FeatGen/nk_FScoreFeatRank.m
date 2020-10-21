function D = nk_FScoreFeatRank(Y, L, N)

% Y     :       Data
% L     :       Target Labels
% L2    :       Nuisance Labels (optional)

ix = unique(L);


if numel(ix) == 2

    indP = L==ix(1); indM = L==ix(2);
    YP = Y(indP,:); YM = Y(indM,:);
    
    % Mean of Positive / Negative Labels
    mP = mean(YP); mM = mean(YM);
    
    % Standard Deviations of Positive / Negative Labels
    sP = std(YP); sM = std(YM);
    vec = [size(Y,2):-1:1]'/size(Y,2);
    if exist('N','var') && ~isempty(N)
        
        ND = zeros(size(Y,2),size(N,2));
        % Compute nuisance-insensitive Fscore(s)
        for i = 1:size(N,2)
            nx = unique(N(:,i));
            indNP = N(:,i)==nx(1); indNM = N(:,i)==nx(2);
            mNP = mean(Y(indNP,:)); mNM = mean(Y(indNM,:));
            sNP = std(Y(indNP,:)); sNM = std(Y(indNM,:));
            tND = (mNP-mNM).^2 ./ (sNP + sNM) ; 
            ND(:,i) = tND';   
        end
    end
    % Compute Fscore and optionally weight down nuisance voxels
    D = (mP-mM).^2 ./ (sP + sM); D = D'; 
    
    if exist('ND','var')
        D = D ./ sum(ND,2);
    end
else
    error('F-Score feature ranking works only for binary problems! Change pre-processing group mode to binary.')
end