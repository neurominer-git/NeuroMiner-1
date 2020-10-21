function [W_ACT, PX] = nk_WAction_config(W_ACT, PX, datadesc, brainmask, defaultsfl, parentstr)

if ~exist('defaultsfl','var') || isempty(defaultsfl),  defaultsfl = 0; end
softflag            = 1;
threshvec           = [25 50 75];
clustfl             = 2; 
exppar              = 1;
act = 0;

if ~defaultsfl
    
    if ~isfield(W_ACT,'softflag'),      W_ACT.softflag = softflag; end
    if ~isfield(W_ACT,'threshvec'),     W_ACT.threshvec = threshvec; end
    if ~isfield(W_ACT,'clustflag'),     W_ACT.clustflag = clustfl; end
    if ~isfield(W_ACT,'exponent'),      W_ACT.exponent = exppar; end
    if ~exist('PX','var'),              PX = []; end
    cluststr            = {'yes','no'};
    weightstr           = {'Weighting','Thresholding'};
    menustr = sprintf('Thresholding or weighting of features using weight vector(s) [ %s ]', weightstr{W_ACT.softflag}); menuact = 1;
    
    switch W_ACT.softflag 
        case 1
            menustr = sprintf('%s|Define exponential multiplier [ %s ]',menustr,nk_ConcatParamstr(W_ACT.exponent)); menuact = [ menuact 4 ];
        case 2
            menustr = sprintf('%s|Define percentile thresholds for feature extraction [ %s ]', menustr, nk_ConcatParamstr(W_ACT.threshvec)); menuact = [ menuact  2];
            if datadesc.type
                menustr = sprintf('%s|Clusterize extracted voxels [ %s ]', menustr, cluststr{W_ACT.clustflag}) ; menuact = [ menuact 3 ];
            end
    end
    
    nk_PrintLogo
    mestr = 'Extract features from rank / weight vector'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 

    act = nk_input(mestr,0,'mq', menustr , menuact);
    
    switch act
        case 1
            if W_ACT.softflag ==1, W_ACT.softflag = 2; else, W_ACT.softflag = 1; end
        case 2
            W_ACT.threshvec = nk_input('Define vector of percentile thresholds',0,'e',W_ACT.threshvec);
        case 3
            W_ACT.clustflag = nk_input('Clusterize voxels based on thresholded weight vector',0,'yes|no',[1,2], W_ACT.clustflag);
            if W_ACT.clustflag, W_ACT.brainmask = brainmask; end
        case 4
            W_ACT.exponent = nk_input('Define exponential multiplier(s)',0,'e', W_ACT.exponent);
          
    end
else 
    W.ACT.softflag = softflag;
    W_ACT.threshvec = threshvec;
    W_ACT.clustflag = clustfl;
    W_ACT.exponent = exppar;
end
PX = nk_AddParam([], [],1, PX,'reset');
switch W_ACT.softflag
    case 1
        PX = nk_AddParam(W_ACT.exponent,'ExpMult', 1, PX);
    case 2
        PX = nk_AddParam(W_ACT.threshvec,'Thresholds', 1, PX);
end
if act, [W_ACT, PX] = nk_WAction_config(W_ACT, PX, datadesc, brainmask, [], parentstr); end
