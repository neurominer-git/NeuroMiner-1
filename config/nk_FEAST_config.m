function [ act, params ] = nk_FEAST_config( params, askFeatNum, parentstr )

NumFeat = 10;
Method= 4;
na_str='?';

if ~exist('askFeatNum','var') || isempty(askFeatNum), askFeatNum = false; end

if isfield(params,'FEAST'), 
    if isfield(params.FEAST,'Method'),  Method = params.FEAST.Method; end
    if isfield(params.FEAST,'NumFeat'), NumFeat= params.FEAST.NumFeat; end
else
    params.FEAST = struct('Method', Method, 'NumFeat', NumFeat);
    params.FEAST = nk_DiscSym_config(params.FEAST,[]);
end

Methods = {'mim','mrmr','cmim','jmi','disr','cife','icap','condred','cmi','mifs'};

if isfield(params.FEAST,'cmd')
    switch params.FEAST.cmd
        case 'discretize'
            disc_str = 'Discretization settings defined';
        case 'symbolize'
            disc_str = 'Symbolization settings defined';
    end
else
    disc_str = na_str;
end

menustr = sprintf('Define discretization settings [ %s ]|', disc_str); menuact = 1;

params.FEAST.MethodStr = Methods{Method};
menustr = [ menustr sprintf('Select FEAST algorithm [ %s ]|', params.FEAST.MethodStr)]; menuact = [menuact 2]; 

if askFeatNum
   FEAST_FeatNum = sprintf('Define number of feature to extract using MRMR [ %g ]|', NumFeat); 
   menustr = [ menustr FEAST_FeatNum ]; menuact = [menuact 3];
end

nk_PrintLogo
mestr = 'FEAST algorithm selection and configuration'; navistr = sprintf('%s\n\t>>> %s',parentstr, mestr); cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
act = nk_input(mestr,0,'mq', menustr, menuact);

switch act
    case 1
        t_act = 1; while t_act > 0, [params.FEAST, t_act] = nk_DiscSym_config( params.FEAST, navistr ); end;
    case 2
        params.FEAST.Method = nk_input('Select FEAST toolbox mutual information algorithm',0,'m', ...
            ['MIM (Mutual Information Maximization; univariate)|', ...
             'MRMR (Maximum Relevance Minimum Redundancy)|', ...
             'CMIM (Conditional Mutual Information Maximization)|', ...
             'JMI (Joint Mutual Information)|', ...
             'DISR (Double Input Symmetrical Relevance)|', ...
             'CIFE (Conditional Infomax Feature Extraction)|', ...
             'ICAP (Interaction Capping)|', ...
             'CONDRED|', ...
             'CMI|', ...
             'MIFS (Mutual Information Feature Selection)|'
             ],1:10, Method);

         switch params.FEAST.Method
             case 10
                 params.FEAST.MethodStr = 'mifs';
                 params.FEAST.MethodParams = nk_input('Define beta (range) for MIFS algorithm',0,'e',1);
             case 11 % BetaGamma not implemented yet
             case 12 % FCBF not implemented yet
         end

         switch params.FEAST.Method
             case {1,2,3,4,5,6,7,8,9}
                 if isfield(params.FEAST,'MethodParams'), params.FEAST = rmfield(params.FEAST,'MethodParams'); end
         end
    case 3
         params.FEAST.NumFeat = nk_input('Number of features to select?',0,'e',NumFeat);
    
end

