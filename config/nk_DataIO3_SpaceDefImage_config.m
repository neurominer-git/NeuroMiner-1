function Thresh = nk_DataIO3_SpaceDefImage_config(Vm, Lm, datatype)

if ~exist('datatype','var') || isempty(datatype), datatype = 'nifti'; end

switch datatype
    case {'nifti','spm'}
        % Read-in space-defining image and check number of unique values
        if ~exist(Vm.fname,'file'), 
            if isempty(fileparts(Vm.fname)), Vm.fname = fullfile(pwd,Vm.fname); end
            if ~exist(Vm.fname,'file'), Thresh = struct('nVml',0,'Vml',0,'Lm',[],'threshop','gt'); return; end
        end
        if ~isstruct(Vm); Vm = spm_vol(Vm); end
        Vmx = spm_read_vols(Vm);
    case 'surf'
        % Read-in space-defining image and check number of unique values
        if isfield(Vm,'fspec') && exist(Vm.fspec,'file'), 
            if isempty(fileparts(Vm.fspec)), Vm.fspec = fullfile(pwd,Vm.fspec); end
            if ~exist(Vm.fspec,'file'), Thresh = struct('nVml',0,'Vml',0,'Lm',[],'threshop','gt'); return; end
            [~,Vmx] = SurfaceReader(Vm.fspec);
        else
            Vmx = Vm.cdata';
        end
end
Vml = unique(Vmx(:)); 
% Ingore 0 label, which should be the background data
if isempty(Vml)
    cprintf('red','\nThis is header-only image. setting threshold to zero.')
    Thresh = struct('nVml',0,'Vml',0,'Lm',[],'threshop','gt');
elseif ~Vml(1), 
    cprintf('red','\nFound a zero label in the space-defining image. Ignoring it.')
    Vml(1)=[]; 
end
nVml = numel(Vml);

% if number of unique values > 1 then we don't have a mask image but a
% label image containing some sort of partitioning information (e.g. an
% atlas image like the AAL atlas)
fprintf('\nDetected %g unique values / labels in the space-defining image', nVml)
if nVml > 1
    ThreshFlag = nk_input('Do you want to threshold the space-defining image?',0,'yes|no',[1,0],1);
    if ThreshFlag
        Vml         = nk_input('Define threshold(s)',0,'e',0); 
        Thresh.threshmode  = nk_input('Define threshold mode',0,'m','Percentile|Absolute',[1,2],2);
        if Thresh.threshmode == 1
            ind0    = Vmx(:)>0;
            Vml     = percentile(Vmx(ind0),Vml);
        end
        nVml        = numel(Vml);
        Lm          = cell(nVml,1);
        threshop    = nk_input('Define threshold operator',0,'m','>|>=|==|<|<=',{'gt','ge','eq','lt','le'},1);
        for i=1:numel(Vml)
            Lm{i}   = sprintf('Binary treshold: %s %g', threshop{1}, Vml(i));
        end
    else
        if ~exist('Lm','var') || isempty(Lm) || numel(Lm)~= nVml
            Lm      = nk_input('Define label (cell array of strings) in ascending Label ID order',0,'e',[],[nVml,1]);
        end
        fprintf('\nSplitting the input space into %g subspaces.',nVml)
        threshop{1} = 'eq';
        
    end
else
    Lm={'Bin Mask'}; threshop{1} = 'gt'; Vml = 0;
end
% Prepare data container
Thresh.nVml     = nVml;
Thresh.Vml      = Vml;
Thresh.threshop = threshop;
Thresh.Lm       = Lm;