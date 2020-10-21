function [param, act] = nk_DiscSym_config(param, parentstr)

dscBinStart     = 0;
dscBinStep      = 0.5;
dscBinStop      = 4;
symMinBin       = 3;
symMaxBin       = 25;
symSeqLength    = 4;
symStdNum       = 3;
cmd             = 'discretize';
methods         = {'discretize', 'symbolize'};

if isfield(param,'DISCRET') && ~isempty(param.DISCRET),
    dscBinStart = param.DISCRET.binstart;
    dscBinStep = param.DISCRET.binsteps;
    dscBinStop = param.DISCRET.binstop;
end

if isfield(param,'SYMBOL') && ~isempty(param.SYMBOL),
    symMinBin = param.SYMBOL.symMinBin;
    symMaxBin = param.SYMBOL.symMaxBin;
    symSeqLength = param.SYMBOL.symSeqLength;
    symStdNum = param.SYMBOL.symStdNum;
end

if isfield(param,'cmd'), cmd = param.cmd; end
switch cmd 
    case 'discretize'
        menustr = [ sprintf('Select binning method [ %s ]|', cmd) ...
                    sprintf('Define binning start value [ %g ]|',dscBinStart) ...
                    sprintf('Define binning stepping [ %g ]|',dscBinStep) ...
                    sprintf('Define binning stop value [ %g ]',dscBinStop) ];
        menuact = 1:4;
        param.DISCRET.binstart = dscBinStart ;
        param.DISCRET.binsteps = dscBinStep;
        param.DISCRET.binstop = dscBinStop;
        if isfield(param,'SYMBOL'), param = rmfield(param,'SYMBOL'); end
    case 'symbolize'  
        menustr = [ sprintf('Select binning method [ %s ]|', cmd) ...
                    sprintf('Define minimum bin count [ %g ]|',symMinBin) ...
                    sprintf('Define maximum bin count [ %g ]|',symMaxBin) ...
                    sprintf('Define sequence length [ %g ]|', symSeqLength) ...
                    sprintf('Define width of data range in number of STDs from mean [ %g ]', symStdNum) ];
        menuact = [1,5:8];
        param.SYMBOL.symMinBin = symMinBin ;
        param.SYMBOL.symMaxBin = symMaxBin ;
        param.SYMBOL.symSeqLength = symSeqLength ;
        param.SYMBOL.symStdNum = symStdNum ;
        if isfield(param,'DISCRET'), param = rmfield(param,'DISCRET'); end
end
param.cmd = cmd;
if ~exist('parentstr','var') || isempty(parentstr), return; end

nk_PrintLogo
mestr = 'Binning algorithm configuration'; navistr = sprintf('%s\n\t>>> %s',parentstr, mestr); cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
act = nk_input(mestr,0,'mq', menustr, menuact);

switch act
            
    case 1
        
        fl = find(strcmp(methods,cmd));
        fl = char(nk_input('Select Binning algorithm',0,'mq', ...
                    ['Unsupervised binning discretization (mean+/-alpha*std)|' ...
                     'Unsupervised entropy-based symbolization'],1:2, fl));
        if fl, param.cmd = methods{fl}; end
        
    case 2
        param.DISCRET.binstart = nk_input('Binning start',0,'r',dscBinStart);
    case 3
        param.DISCRET.binsteps = nk_input('Binning step width (alpha*std)',0,'r',dscBinStep);
    case 4
        param.DISCRET.binstop = nk_input('Binning stop',0,'r',dscBinStop);      
    case 5
        param.SYMBOL.symMinBin = nk_input('Minimum bin count',0,'i',symMinBin);
    case 6
        param.SYMBOL.symMaxBin = nk_input('Maximum bin count',0,'i',symMaxBin);
    case 7
        param.SYMBOL.symSeqLength = nk_input('Sequence length',0,'i',symSeqLength);
    case 8
        param.SYMBOL.symStdNum = nk_input('Width of data range (# of STDs from mean)',0,'r',symStdNum);
                         
end

