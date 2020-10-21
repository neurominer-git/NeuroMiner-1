function PP = nk_AddParam(p, desc, t, PP, actstr)

if ~exist('actstr','var'), actstr = []; end
if strcmp(actstr,'reset')
    Px = struct('Params',[],'Params_desc',[],'Typ',[]); PP.Px = Px; PP.opt = []; return
elseif ~exist('PP','var') || isempty(PP) || strcmp(actstr,'replace');
    Px = struct('Params',[],'Params_desc',[],'Typ',[]);
elseif isfield(PP,'Px')
    Px = PP.Px;    
end

nP=1; multfl = 0;
if iscell(p) 
    if numel(p) > 1
        nP = numel(p); multfl = 1; 
    else
        p = p;
    end
end

for j=1:nP
    if multfl, 
        tp = p{j}; tdesc = desc{j}; tt = t{j};
    else
        tp = p; tdesc = desc; tt = t;
    end
    % Check whether parameter already exists
    fl = 0;
    for i=1:numel(Px)
        if strcmp(Px(i).Params_desc, tdesc), fl = 1; break; end
    end
    if fl
        if isempty(tp)
            % Delete
            Px(i) = [];
        else
            % Overwrite
            Px(i).Params         = tp;
            Px(i).Params_desc    = tdesc;
            Px(i).Typ            = tt;
        end
    else
        ll = numel(Px); if ~isempty(Px(ll).Params), ll=ll+1; end 
        % Add to end
        Px(ll).Params         = tp;
        Px(ll).Params_desc    = tdesc;
        Px(ll).Typ            = tt;
    end
end
cPx = struct2cell(Px); cParams = squeeze(cPx(1,1,:));
PP.opt = allcomb2(cParams,'matlab'); 
PP.Px = Px;