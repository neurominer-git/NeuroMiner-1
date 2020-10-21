function imrelief = nk_IMRelief_config(imrelief, classhdr, setupfl, defaultsfl)
% function param = nk_IMRelief_config(param, setupfl)
%
% This function sets up the parameters for the IMRelief feature 
% selection algorithm presented by Yijun Sun and colleagues:
%
% IMRelief_Sigmoid_FastImple: fast implemetation of IMRelief_Sigmoid
% --------------------------------------------------------------------------
% INPUT:
%       patterns:  training data: [x1,x2,...xn] Each column is an observation
%        targets:  class label = {1,2}
%           Para:  parameters. 
%  Para.distance:  distance metric (default: block distance)
%     Para.sigma:  kernel width
%    Para.lambda:  regulariztion parameter
%      Para.plot:  1: plot of the learning process; 0: do not plot
%   Para.maxiter:  max iterations to minimize theta
%
% OUTPUT:
%         Weight:  weight of features
% -------------------------------------------------------------------------
% by Yijun Sun @University of Florida 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) N. Koutsouleris 04/2015

% Default Settings
dist = 1; diststr = 'block';
sigma = 2;
lambda = 0.1;
plt = 0;
maxiter = 10;
CV2flag = 1;
if ~exist('classhdr','var'), classhdr='IMRelief as classification algorithm'; end
if ~exist('setupfl','var'), setupfl=0; end
if ~exist('defaultsfl','var'), defaultsfl=0; end
if ~setupfl && exist('imrelief','var') && ~isempty(imrelief)
    if isfield(imrelief,'distance')
        diststr = imrelief.distance;
        switch diststr
            case 'block'
                dist=1';
            case 'euclidean'
                dist=2;
        end
    end
    if isfield(imrelief,'sigma')
        sigma = imrelief.sigma;
    end
    if isfield(imrelief,'lambda')
        lambda = imrelief.lambda;
    end
    if isfield(imrelief,'maxiter')
        maxiter = imrelief.maxiter;
    end
%     if isfield(imrelief,'CV2flag')
%         CV2flag = imrelief.CV2flag + 1;
%     end
end
if ~defaultsfl
    % -------------------------------------------------------------------------
    nk_PrintLogo
    lambdastr = nk_ConcatParamstr(lambda);
    sigmastr = nk_ConcatParamstr(sigma);
    
    menuact = ['Distance measure [ ' diststr ' ]|' ...
    'Kernel parameter ''Sigma'' defining resolution of local learning [ ' sigmastr ' ]|' ...
    'Regularization parameter ''Lambda'' [ ' lambdastr ' ]|' ...
    'Maximum iterations [ ' num2str(maxiter) ' ]'];
    
    menusel = 1:4;

%     if numel(sigma) > 1 || numel(lambda) > 1
%        menuact = [menuact ...
%            '|Optimization at single or across all CV1 partitions' ];
%        menusel = [menusel 5];
%     end
    
    act = nk_input(['IMRelief: Parameter setup for ' classhdr ],0, 'mq', menuact, menusel);
    
    switch act
        case 1
            dist = nk_input('Distance measure',0,'mq', ...
                ['Block distance|' ...
                'Euclidean distance'], [1, 2], dist);
        case 2
            sigma = nk_input('Sigma',0,'e',sigma);
        case 3
            lambda = nk_input('Lambda',0,'e',lambda);
        case 4
            maxiter = nk_input('Maximum iterations',0,'w1',maxiter);
        case 5
            CV2flag = nk_input('Optimization mode',0,'mq', ...
                'Individually|Across all CV1 partitions',1:2, CV2flag);
            if CV2flag
                CV2flag = CV2flag - 1;
                imrelief.CV2flag = CV2flag;
            end
    end
else
    act = 0;
end
switch dist
    case 1
        imrelief.distance = 'block';
    case 2
        imrelief.distance = 'euclidean';
end
imrelief.desc = classhdr;
imrelief.sigma = sigma;
imrelief.lambda = lambda;
imrelief.maxiter = maxiter;
imrelief.plotfigure = plt;
if act, imrelief = nk_IMRelief_config(imrelief, classhdr); end

return