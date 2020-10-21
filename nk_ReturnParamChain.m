function PX = nk_ReturnParamChain(PO, optimonly)

if ~exist('optimonly','var') || isempty(optimonly), optimonly = false; end
nP = 1; if iscell(PO), nP = numel(PO); end
PX = struct('Params',[],'Params_desc',[],'Typ',[],'cmd',[],'steps',[],'stepparams',[],'opt',[]);

for n=1:nP
    
    Params = []; Params_desc = []; Typ = []; cmd = ''; steps = []; ll=1;
    if iscell(PO), PREPROC = PO{n}; else PREPROC = PO; end
    if nP>1, Pstr = sprintf('_M%g',n); end
    if isfield(PREPROC,'PX') && ~isempty(PREPROC.PX) && isfield(PREPROC.PX,'Px')
        [Params, Params_desc, Typ] = retrieve_from_struct(PREPROC.PX.Px, optimonly);
        if isfield(PREPROC,'cmd'), cmd = PREPROC.cmd; else cmd=[]; end
        steps = 1;
        f = fieldnames(PREPROC);
        stepparams = getfield(PREPROC,f{1});
    else
        % First check processing sequence for optimisable parameters
        if isfield(PREPROC,'ACTPARAM')
            % loop from back to beginning !!!
            for i=numel(PREPROC.ACTPARAM):-1:1 
               if isfield(PREPROC.ACTPARAM{i},'PX') && ~isempty(PREPROC.ACTPARAM{i}.PX) && isfield(PREPROC.ACTPARAM{i}.PX,'Px')
                   [cParams, cParams_desc, cTyp] = retrieve_from_struct(PREPROC.ACTPARAM{i}.PX.Px, optimonly);
                   if isempty(cParams), continue; end
                   Params = [Params cParams];
                   Params_desc = [Params_desc cParams_desc];
                   Typ = [Typ cTyp];
                   cmd = char(cmd, repmat(PREPROC.ACTPARAM{i}.cmd,numel(cParams_desc),1));
                   steps = [steps repmat(i,1,numel(cParams_desc))];
                   %f = fieldnames(PREPROC.ACTPARAM{i});
                   stepparams{ll} = PREPROC.ACTPARAM{i};
                   ll = ll+1;
               end
            end
        end
        % Finally, check for pre-chain parameter arrays
        if isfield(PREPROC,'SPATIAL') && isfield(PREPROC.SPATIAL,'PX') && ~isempty(PREPROC.SPATIAL.PX.opt)
             [cParams, cParams_desc, cTyp] = retrieve_from_struct(PREPROC.SPATIAL.PX.Px, optimonly);
             if ~isempty(cParams), 
                 Params = [Params cParams];
                 Params_desc = [Params_desc cParams_desc];
                 Typ = [Typ cTyp];
                 cmd = strvcat(cmd, 'spatialfilter');
                 steps = [steps ll];
                 stepparams{ll} = PREPROC.SPATIAL;
             end
        end
    end

    if ~isempty(Params)
        if iscell(PO)
            PX(n).Params = Params;
            PX(n).Params_desc = Params_desc;
            PX(n).Typ = Typ;
            if ~isempty(cmd), 
                PX(n).cmd = cellstr(cmd); 
                 if strcmp(PX(n).cmd{1},''), PX(n).cmd(1)=[]; end
            end
            PX(n).steps = steps;
            PX(n).stepparams = stepparams;
            PX(n).opt = allcomb2(PX(n).Params,'matlab'); 
        else
            PX.Params = Params;
            PX.Params_desc = Params_desc;
            PX.Typ = Typ;
            if ~isempty(cmd), 
                PX.cmd = cellstr(cmd);
                if strcmp(PX.cmd{1},''), PX.cmd(1)=[]; end
            end
            PX.steps = steps;
            PX.stepparams = stepparams;
            PX.opt = allcomb2(PX(n).Params,'matlab'); 
        end
    else
        if ~iscell(PO); PX = []; end
    end
end

function [Params, Params_desc, Typ] = retrieve_from_struct(Px, optimonly)

   Px = struct2cell(Px);
   if optimonly
       indparam = cellfun(@numel,squeeze(Px(1,1,:))) > 1;
   else
       indparam = 1:size(Px,3);
   end
   Params = squeeze(Px(1,1,indparam))';
   Params_desc = squeeze(Px(2,1,indparam))';
   Typ = squeeze(Px(3,1,indparam))';

