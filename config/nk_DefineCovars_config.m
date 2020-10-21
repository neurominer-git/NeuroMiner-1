function [G, Gnames] = nk_DefineCovars_config(n_subjects_all, covars)

if exist('covars','var') && ~isempty(covars),
    ncov = size(covars,2);
else
    ncov = Inf;
end

G = nk_input('Covariates',0,'r',[],[n_subjects_all ncov]);

if nargout == 2
    Gnames =  cell(size(G,2),1);
    descflag = nk_input('Define covariate names',0,'no|one-by-one|vector',[0,1,2],0);
    switch descflag
        case {0,1}

            covstr = sprintf('\nFirst row of covariate matrix:');
            covstr = [covstr sprintf('\t%g',G(1,:))]; 
            disp(covstr)

            for i=1:size(G,2)
                if descflag
                    Gnames{i} = nk_input(['Provide name of covariate ' num2str(i)],0,'s');
                else
                    Gnames{i} = ['covar' num2str(i)];
                end
            end
        case 2
            Gnames = nk_input('Covariate name vector',0,'e',[],size(G,2));
            if ~iscell(Gnames) 
                if isnumeric(Gnames)
                    error('Covariate names cannot be numeric but have to be a cell array of strings')
                else
                    Gnames = cellstr(Gnames);
                end
            end
    end
else
    Gnames = [];
end

return