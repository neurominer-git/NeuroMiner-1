function PREPROC = nk_SetParamChain(paramfl, curclass, PREPROC)

if isfield(paramfl,'PREPROC') && isfield(paramfl,'PXfull') && ~isempty(paramfl.PXfull)
    PREPROC = paramfl.PREPROC; %cnt = numel(paramfl.PREPROC.ACTPARAM); curstep = 0;
    for y=1:size(paramfl.PXfull.cmd,1)
        % Set ACTPARAM
        if ~exist('PX','var') || paramfl.PXfull.steps(y) ~= curstep
            PX = nk_AddParam([], [], [], [], 'reset'); 
            curstep = paramfl.PXfull.steps(y);
        end
         % Finally, check for pre-chain parameter arrays
        if strcmp(paramfl.PXfull.cmd(y,:),'spatialfilter')
             if ~exist('PX','var'),  PX = nk_AddParam([], [], [], [], 'reset');  end
             PREPROC.SPATIAL.PX =  nk_AddParam( unique(paramfl.PXunique{curclass}(:,y)), paramfl.PXfull.Params_desc{y}, 1, PX, 'replace');
        else
            % Loop through every parameter stored in PREPROC.ACTPARAM.PX.Px
            % and check whether it has been optimized as defined in
            % paramfl.PXfull. If so, write optimized parameter into Px,
            % otherwise take the parameter defined in PREPROC
            for z=1:numel(PREPROC.ACTPARAM{paramfl.PXfull.steps(y)}.PX.Px)
                if strcmp(paramfl.PXfull.Params_desc{y}, PREPROC.ACTPARAM{paramfl.PXfull.steps(y)}.PX.Px(z).Params_desc)
                    PX = nk_AddParam( unique(paramfl.PXunique{curclass}(:,y))', paramfl.PXfull.Params_desc{y}, 1, PX );  
                else
                    PX = nk_AddParam( PREPROC.ACTPARAM{paramfl.PXfull.steps(y)}.PX.Px(z).Params, PREPROC.ACTPARAM{paramfl.PXfull.steps(y)}.PX.Px(z).Params_desc, 1, PX );  
                end
            end
            PREPROC.ACTPARAM{paramfl.PXfull.steps(y)}.PX = PX; 
        end
    end
end