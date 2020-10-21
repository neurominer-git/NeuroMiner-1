function [tY, Fx] = nk_ExtractFeatures(Y, F, Ytest, indF)

if iscell(Y)
    nvar    = size(Y,2);
    tY      = cell(1,nvar);
    for v=1:nvar
        if isempty(Ytest)
            if isempty(F)
                tY{v} = Y{v};
            else
                D = size(Y{v},2);
                Fx = prep_fmask(D, F{v}, indF);
                tY{v} = Y{v}(:,Fx);
            end
        else % Concatenate Training and Test data
            if isempty(F)
                tY{v} = [Y{v}; Ytest{v}];
            else
                Fx = prep_fmask(D, F{v}, indF);
                tY{v} = [Y{v}(:,Fx); Ytest{v}(:,Fx)];
            end
        end
    end
    if nvar == 1, tY = tY{1}; end
else
    if isempty(Ytest)
        
        if isempty(F)
            tY = Y;
        else
            D = size(Y,2);
            Fx = prep_fmask(D, F, indF);
            tY = Y(:,Fx);
        end
    else
        if isempty(F)
            tY = [Y; Ytest];
        else
            D = size(Y,2);
            Fx = prep_fmask(D, F, indF);
            tY = [Y(:,Fx); Ytest(:,Fx)]; 
        end
    end
end

return

function Fx = prep_fmask(D, F, indF)

if size(F,2) == 1 && sum(any(F,2)) == numel(F)
    
    Fx = F;
    
else

    Fx = logical(F(1:D,indF));

end

return