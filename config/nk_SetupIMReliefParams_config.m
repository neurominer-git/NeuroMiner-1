function param = nk_SetupIMReliefParams_config(res, param, setupfl)

if ~exist('setupfl','var')
    if isfield(param,'imrelief'), setupfl = false; else setupfl = true; end
end

nclass = numel(unique(res.label));
if  nclass > 2
    fprintf('\nNeuroMiner detected that you have specified %g binary classifiers in your NM structure.', nclass )
    binflg = nk_input('Setup separate IMRelief parameters for each of these binary comparisons',0,'yes|no',[1,0],1);
    if ~binflg, 
        nclass = 1;
        if isfield(param,'imrelief') && iscell(param.imrelief) && numel(param.imrelief) > 1
            imrelief = param.imrelief{1};
            param.imrelief = [];
            param.imrelief{1} = imrelief;
        end
    end
else
    binflg = true; nclass  = 1;
end

for curclass = 1 : nclass
    if ~isfield(param,'imrelief')
        imrelief = [];
        classhdr = sprintf('for all binary classifiers');
    else
        if iscell(param.imrelief)
            if numel(param.imrelief) < curclass
                imrelief = param.imrelief{end};
            else
                imrelief = param.imrelief{curclass};
            end
        else
            imrelief = param.imrelief;
            param.imrelief = [];
        end

        if binflg, 
            classhdr = sprintf('binary classifier %s', res.cv.class{1,1}{curclass}.groupdesc);
        else
            classhdr = sprintf('for all binary classifiers'); 
        end
    end
    param.imrelief{curclass} = nk_IMRelief_config(imrelief, classhdr, setupfl, 0);
end
    
end