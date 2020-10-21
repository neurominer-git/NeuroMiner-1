function nk_ThreshWriteImg(brainmask, P, resl)

if ~exist('P','var') || ~exist(deblank(P(1,:)),'file')
    P = spm_select(Inf,'image','Select images for thresholding');
end
if ~exist('brainmask','var') || ~exist(brainmask,'file')
    brainmask = spm_select(1,'image','Brainmask image');
    brainmask = regexprep(brainmask,',1','');
end
if ~exist('resl','var') || ~exist(resl,'file')
    reslflag = nk_input('Reslice output',0,'yes|no',[1,0],1);
    if reslflag
        resl = spm_select(1,'image','Space-defining image');
        resl= regexprep(resl,',1','');
    end
else
    reslflag = true;
end

if reslflag 
    prec = 'uint8|int16|int32|float32|float64|int8|uint16|uint32';
    types   = [    2      4      8   16   64   256    512    768];
    datatype = nk_input('Select data type',0,'m',prec,types);
end
nP = size(P,1);

thresh = config_threshold;

sclflag = nk_input('Scale image?',0,'yes|no',[1,0],0);
if sclflag
   sclminmax = nk_input('Minimum / Maximum values',0,'e');
   sclsuff = '_scl';
else
    sclsuff = '';
end

for i=1:nP
   Px = deblank(P(i,:));
   Px = regexprep(Px,',1','');
   [dims, indvol, Yi] = nk_ReadMaskIndVol(brainmask,Px);
   if sclflag
       fprintf('\nScaling image to range [ %g <-> %g ]',sclminmax(1),sclminmax(2));
       Yi = nk_ScaleData(Yi,sclminmax(1),sclminmax(2));
   end
   [p,n,e] = fileparts(Px);
   [T, Thresh] = nk_Threshold(Yi, thresh);
   if ~isempty(Thresh)
       MaxT = max(T); MinT = min(T);
       Vec =[ Thresh, MaxT, MinT];
       fprintf('\nThresholds:')
       fprintf('\t%g',Vec);
       filnam = [n sclsuff '_thresh'];
   else
       filnam = n;
   end
   nk_WriteVol(T, filnam, 2, brainmask);
   if reslflag
       if ~isempty(Thresh)
           filnam2 = [n sclsuff '_thresh_resl'  e];
       else
           filnam2 = [n sclsuff '_resl'  e];
       end
       Inp = char(resl,[filnam e]);
       spm_imcalc_ui(Inp,filnam2,'i2',{[],[],datatype,1});
   end
end
fprintf('\nDone.\n')

end


function thresh = config_threshold(threshtype, threshval, threshlogop)

def = [];
if exist('threshtype','var') && ~isempty(threshtype)
    def = threshtype;
end

thresh.type = nk_input(['Threshold type'],0,'m', ...
                            ['Percentile|' ...
                            'Absolute Value|' ...
                            'None'],1:3,def);
                        
def = [];
if exist('thresval','var') && ~isempty(threshval)
    def = threshval;
end

switch thresh.type
    case 1
        thresh.val = nk_input(['Define percentile(s) for threshold'],0,'e',def);
    case 2
        thresh.val = nk_input(['Define absolute value(s) for threshold'],0,'e',def);
    case 3
        thresh.val = [];
end


if thresh.type == 1 || thresh.type == 2
    def = [];
    if exist('threslogop','var') && ~isempty(threshlogop) && ~threshlogop
        def = threshlogop;
    end
    switch length(thresh.val)
        case 1
            thresh.logop = nk_input('Define logical operation for thresholding', 0, 'm', ...
                ['>|' ...
                '>=|' ...
                '<|' ...
                '<='],1:4, def);
        case 2
            thresh.logop = nk_input('Define logical operation for thresholding', 0, 'm', ...
                ['>< (less than ' num2str(thresh.val(1)) ' OR greater than ' num2str(thresh.val(2)) ')|' ...
                '>=< (less or equal than ' num2str(thresh.val(1)) ' OR greater or equal than ' num2str(thresh.val(2)) ')|' ...
                '<> (greater than ' num2str(thresh.val(1)) ' and less than ' num2str(thresh.val(2)) ')|' ...
                '<=> (greater or equal than ' num2str(thresh.val(1)) ' and less or equal than ' num2str(thresh.val(2)) ')'],[7,8,5,6], def);
    end
else
    thresh.logop = 0;
end

end