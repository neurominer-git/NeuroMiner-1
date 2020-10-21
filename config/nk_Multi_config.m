function [MULTI, act] = nk_Multi_config(MULTI, defaultfl, parentstr)
% function res = nk_MULTI_config(res)
%
% Setup parameters for multi-group classification
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) N. Koutsouleris 04/2015
global NM

if ~exist('defaultfl','var') || isempty(defaultfl), defaultfl = 0; end

% Defaults:
% =========
if numel(unique(NM.label)) < 3
    MULTI.flag = 0;
    MULTI.train = 0;
    return; 
end
multiflag     = 1;
hardcoded     = 1;
multitrain    = 1;
method        = 2;
coding        = 1;
decoding      = 1;
BinBind       = 0;
decisiontype  = 1;
act           = 0;
if ~defaultfl
    nk_PrintLogo
    % Take over previous values if available
    if ~isempty(MULTI)
        if isfield(MULTI,'flag'),       multiflag = MULTI.flag; end
        if isfield(MULTI,'hardcoded'),  hardcoded = MULTI.hardcoded; end
        if isfield(MULTI,'train'),      multitrain = MULTI.train; end
        if isfield(MULTI,'method'),     method = MULTI.method; end    
        if isfield(MULTI,'coding'),     coding= MULTI.coding; end
        if isfield(MULTI,'decoding'),   decoding = MULTI.decoding; end
        if isfield(MULTI,'BinBind'),    BinBind = MULTI.BinBind; end
    end
    
    if ~multiflag,  
        multiflagstr = 'no';    
        menustr = sprintf('Train multi-class predictor [ %s ]', multiflagstr) ;
        menuact = 1;
        multiflagdef = 2;
    else
        multiflagdef = 1;
        multiflagstr = 'yes';
        if BinBind, binbinddef = 1; else binbinddef = 2; end
        if ~multitrain, 
            multitrainstr =  'no'; 
            multitraindef = 2;
        else
            multitrainstr = 'yes';
            multitraindef = 1;
        end
        switch method
            case 1
                multimethodstr = 'One-vs-One-Max-Wins ';
                switch MULTI.decisiontype
                    case 1
                        decodestr = 'Sum';
                    case 2
                        decodestr  = 'Mean';
                    case 3
                        decodestr = 'Product';
                    case 4
                        decodestr = 'Majority';
                    case 5
                        decodestr = 'Median';
                end
                
            case 2
                multimethodstr = 'Error Correcting Output Codes ';
                switch MULTI.decoding
                    case 1
                        decodestr = 'Hamming distance';
                    case 2
                        decodestr = 'Euclidean distance';
                    case 3
                        decodestr = 'Laplacian decoding';
                    case 4
                        decodestr = 'Attenuated euclidean distance';
                    case 5
                        decodestr = 'Linear loss-based decoding';
                end
            case 3
                multimethodstr = 'Directed Acyclic Graph '; decodestr = 'DAG';
        end
        multimethodstr = [multimethodstr '( ' decodestr ' )'] ;
        menustr = [ sprintf('Train multi-class predictor [ %s ]|', multiflagstr) ...
                sprintf('Optimize NM training process for multi-group classification performance [ %s ]|', multitrainstr) ...
                sprintf('Specify multi-class decision mechanism [ %s ]', multimethodstr) ];
        menuact = 1:3;
    end        
    
    mestr = 'Multi-class prediction parameters'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\n\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq',menustr, menuact);
    
    switch act
        case 1
            multiflag = nk_input('Train multi-class predictor',0,'yes|no',[1,0],multiflagdef);
        case 2
            multitrain = nk_input('Optimize NM training process for multi-group classification performance',0,'yes|no',[1,0],multitraindef);
            if ~multitrain, BinBind = nk_input('Bind multi-class predictor to binary classifiers'' optima',0,'yes|no',[1,0], binbinddef); else BinBind = 0; end
        case 3
             if isfield(NM.TrainParam,'RAND') && ...
                isfield(NM.TrainParam.RAND,'Decompose') && ...
                    NM.TrainParam.RAND.Decompose == 2

                method = nk_input('Multi-class decision method',0,'m',...
                    ['Simple Pairwise Decoding (Maximum-Wins method)|' ...
                    'Error-correcting output codes'],1:2,method);
            else

                method = nk_input('Multi-class decision method',0,'m',...
                    ['Simple Pairwise Decoding (Maximum-Wins method)|' ...
                    'Error-correcting output codes|' ...
                    'Directed Acyclic Graph'],1:3,method);
             end
             switch method
                case 1
                    decisiontype = nk_input('Multi-class decision type',0,'m', ...
                        ['Sum of decision values|' ...
                         'Mean of decision values|' ...
                         'Product of decision values|' ...
                         'Majority voriting|' ...
                         'Median of decision values'],1:5,decisiontype);
                case 2
                    coding = 1; % Pairwise
                    decoding = nk_input('Decoding method',0,'m', ...
                        ['Hamming distance|' ...
                        'Euclidean distance|' ...
                        'Laplacian decoding|' ...
                        'Attenuated euclidean distance|' ...
                        'Linear loss-based decoding'], 1:5,decoding);
            end
    end
end
MULTI.flag           = multiflag;
MULTI.method         = method;
MULTI.train          = multitrain;
MULTI.hardcoded      = hardcoded;
MULTI.coding         = coding;
MULTI.decoding       = decoding;
MULTI.BinBind        = BinBind;
MULTI.decisiontype   = decisiontype;
