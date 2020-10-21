function [varind, varstr] = nk_SelectVariateIndex(dat, unifl, varind, askfl, showvar)
% ----------------------------------------------------------------------------------
% FORMAT [varind, varstr] = nk_SelectVariateIndex(dat, unifl, varind, askfl, showvar)
% ----------------------------------------------------------------------------------
% Select modality to work on depending on dat.Y
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 09 / 2017

nvar = length(dat.Y); 
if ~exist('unifl','var') ||  isempty(unifl),    unifl = 1;  end
if ~exist('varind','var') || isempty(varind),   varind = 1; end
if ~exist('askfl','var') || isempty(askfl),     askfl = 1;  end
if ~exist('showvar','var') || isempty(showvar), showvar = 1:nvar;  end

if iscell(dat.Y),
           
    if isfield(dat,'datadescriptor')
        
        nk_PrintLogo
        [nanCL, nCL] = nm_nancount(dat.label);
        fprintf('\n'); cprintf('black*','NM framework: %s ==> ', dat.modeflag);
        switch dat.modeflag
            case 'regression'
               fprintf('%g cases, median (SD): %g (%g)', nanCL, nm_nanmedian(dat.label), nm_nanstd(dat.label));
            case 'classification'
               fprintf('%g cases in %g samples ( %s )', nanCL, numel(dat.n_subjects), strjoin(dat.groupnames,', '));               
        end
        if nCL>0, fprintf('; %g cases without labels', nCL); end
        fprintf('\n')
        cprintf('blue','\n*****************************')
        cprintf('blue','\n*****    MODALITIES     *****')  
        cprintf('blue','\n*****************************')
        fprintf('\n')
        menuact = []; ll = 1;
        for i=1:nvar
            if isempty(find(showvar == i)), continue; end
            menuact = [menuact i]; badcases = []; badfeats = [];
            [~, ~, ~, NanStats] = nk_ManageNanCases(dat.Y{i}, [], [], 'inform');
            if sum(NanStats.Cases)>0,
                badcases = sprintf('%1.1f%% cases with NaNs', sum(NanStats.Cases>0)*100/numel(NanStats.Cases));
                if any(NanStats.Cases100), badcases = sprintf('%s | %g (%1.1f%%) = 100%% NaNs', badcases, sum(NanStats.Cases100),NanStats.Cases100Perc); end
                if any(NanStats.Cases50), badcases = sprintf('%s | %g (%1.1f%%) = [50%%-99%%] NaNs', badcases, sum(NanStats.Cases50),NanStats.Cases50Perc); end
                if any(NanStats.Cases25), badcases = sprintf('%s | %g (%1.1f%%) = [25%%-50%%] NaNs', badcases, sum(NanStats.Cases25),NanStats.Cases25Perc); end
            end
            if sum(NanStats.Feats)>0,
                badfeats = sprintf('%g%% features with NaNs', sum(NanStats.Feats>0)*100/numel(NanStats.Feats));
                if any(NanStats.Feats100), badfeats = sprintf('%s | %g (%1.1f%%) = 100%% NaNs', badfeats, sum(NanStats.Feats100),NanStats.Feats100Perc); end
                if any(NanStats.Feats50), badfeats = sprintf('%s | %g (%1.1f%%) = [50%%-99%%] NaNs', badfeats, sum(NanStats.Feats50),NanStats.Feats50Perc); end
                if any(NanStats.Feats25), badfeats = sprintf('%s | %g (%1.1f%%) = [25%%-50%%] NaNs', badfeats, sum(NanStats.Feats25),NanStats.Feats25Perc); end
            end
            fprintf('\n[ %3g ]',i);cprintf('black*',' %s', dat.datadescriptor{i}.desc);
            if nvar < 10
                fprintf('\n');
                switch dat.datadescriptor{i}.input_settings.datasource
                    case 'matrix'
                        switch dat.datadescriptor{i}.input_settings.groupmode
                            case 1
                                fprintf('\t\t* Source: MATLAB workspace\n')
                            otherwise
                                switch dat.datadescriptor{i}.input_settings.groupmode
                                    case 4
                                        sheet_str = sprintf(' [ Table: %s ]',dat.datadescriptor{i}.input_settings.sheet);
                                    otherwise
                                        sheet_str = '';
                                        
                                end
                                fprintf('\t\t* Source: %s: %s%s\n',dat.datadescriptor{i}.input_settings.M_filestr, dat.datadescriptor{i}.input_settings.M_edit, sheet_str);     
                        end
                    case {'spm','nifti'}
                        switch dat.datadescriptor{i}.input_settings.datasource
                            case 'nifti'
                                fprintf('\t\t* Source: NIFTI/ANALYZE images\n')
                            case 'spm'
                                fprintf('\t\t* Source: SPM design matrix\n')
                        end
                        if ~isempty(dat.brainmask{i}),
                            fprintf('\t\t* Space-defining image: %s\n', dat.brainmask{i});
                        end
                        if isfield(dat.datadescriptor{i},'globnorm') 
                           fprintf('\t\t* Global normalization factor: %g\n', dat.datadescriptor{i}.globnorm)
                        end
                        if isfield(dat.datadescriptor{i},'globscale')
                            fprintf('\t\t* Global scaling option: ')
                            switch dat.datadescriptor{i}.globscale
                                case 1
                                    fprintf('not applied.\n')
                                case 2
                                    fprintf('user-defined global scaling applied.\n')
                                    if dat.datadescriptor{i}.input_settings.gopt == 2
                                       fprintf('\t\t* Globals file: %s\n',dat.datadescriptor{i}.input_settings.globvar_edit)
                                    end
                                case 3
                                    fprintf('automatic global scaling applied.\n')
                            end
                        end
                        
                    case 'surf'
                        fprintf('\t\t* Source: Surface files\n')
                end
                fprintf('\t\t* Dimensionality: %g features %s\n', size(dat.Y{i},2))
                if ~isempty(badcases), fprintf('\n\t\t'); cprintf([1,0.5,0],'* NaNs per cases: %s ', badcases); end
                if ~isempty(badfeats); fprintf('\n\t\t'); cprintf([1,0.5,0],'* NaNs per features: %s ', badfeats); end
            end
            ll=ll+1;
        end
    end
    
    if askfl && nvar > 1 && nargout > 0
    
        if unifl 
            varind = nk_input('Specify a modality to work on',0,'i', 1);
        else
            varind = nk_input('Specify one or more modalities to work on',0, 'i');
        end
        if varind <= numel(menuact)
            varind = menuact(varind); 
        end
    end
   
end

varstr = ['_var' num2str(varind)];

return