function [ act, params ] = nk_MRMR_config( params, askFeatNum, parentstr )

mRMRscheme = 1;
NumFeat = 10;
na_str = '?';

MRMRschemes = {'MID','MIQ'};
if ~exist('askFeatNum','var') || isempty(askFeatNum), askFeatNum = false; end

if isfield(params,'MRMR')
    if isfield(params.MRMR,'mRMRscheme'), mRMRscheme = params.MRMR.mRMRscheme; end
    if isfield(params.MRMR,'NumFeat'), NumFeat= params.MRMR.NumFeat; end
else
    params.MRMR = struct('mRMRscheme', mRMRscheme, 'NumFeat', NumFeat);
    params.MRMR = nk_DiscSym_config(params.MRMR,[]);
end

if isfield(params.MRMR,'cmd')
    switch params.MRMR.cmd
        case 'discretize'
            disc_str = 'Discretization settings defined';
        case 'symbolize'
            disc_str = 'Symbolization settings defined';
    end
else
    disc_str = na_str;
end

menustr = sprintf('Define discretization settings [ %s ]|', disc_str); menuact = 1;

scheme_str = MRMRschemes{mRMRscheme}; 
menustr = [ menustr sprintf('Define MRMR scheme [ %s ]|', scheme_str)]; menuact = [menuact 2]; 

if askFeatNum
   MRMR_FeatNum = sprintf('Define number of feature to extract using MRMR [ %g ]|', NumFeat); 
   menustr = [ menustr MRMR_FeatNum ]; menuact = [menuact 3];
end

if ~exist('parentstr','var') || isempty(parentstr), return; end

nk_PrintLogo
mestr = 'Maximum Relevance Minimum Redundancy (MRMR) algorithm configuration'; navistr = sprintf('%s\n\t>>> %s',parentstr, mestr); cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
act = nk_input(mestr,0,'mq', menustr, menuact);

switch act 
    case 1
        t_act = 1; while t_act > 0, [params.MRMR, t_act] = nk_DiscSym_config( params.MRMR, navistr ); end;
    case 2
        params.MRMR.mRMRschteme = nk_input('mRMR scheme',0,'MID|MIQ',[1,2], mRMRscheme);
    case 3
        params.MRMR.NumFeat = nk_input('Number of features to select?',0,'e',NumFeat);
end

