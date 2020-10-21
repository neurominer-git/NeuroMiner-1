function param = nk_RGS_config(param, setupfl, defaultsfl)
% function param = nk_RGS_config(param, setupfl)
%
% Setup parameters for the RGS algorithm (kNN-regression)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) N. Koutsouleris 12/2010

nperms = 50;
nMinFeatsPerc = 10;
nMaxFeatsPerc = 90;

if ~exist('setupfl','var'), setupfl=0; end
if ~exist('defaultsfl','var'),defaultsfl=0; end

if ~defaultsfl % use interactive setup
    if ~setupfl
        if isfield(param,'RSS')
            if isfield(param.RSS,'nperms'),        nperms        = param.RSS.nperms; end
            if isfield(param.RSS,'nMinFeatsPerc'), nMinFeatsPerc = param.RSS.nMinFeatsPerc; end
            if isfield(param.RSS,'nMaxFeatsPerc'), nMaxFeatsPerc = param.RSS.nMaxFeatsPerc; end
        end
    end

    % -------------------------------------------------------------------------
    nk_PrintLogo

    act = nk_input('Ransom Subspace Sampling: Parameter Setup',0, 'm', ...
            ['# Permutations [' num2str(nperms) ']|' ...
            'Minimum percentage of features per subspace [' num2str(nMinFeatsPerc) ']|' ...
            'Maximum percentage of features per subspace [' num2str(nMaxFeatsPerc) ']|' ...
            '<< Back'],1:4);

    switch act
        case 1
            nperms = nk_input('# Permutations',0,'e', nperms);
        case 2
            nMinFeatsPerc = nk_input('Minimum percentage of features per subspace',0,'e',nMinFeatsPerc);
        case 3
            nMaxFeatsPerc = nk_input('Maximum percentage of features per subspace',0,'e',nMaxFeatsPerc);
    end
else
    act = 4;
end
param.RSS.nperms        = nperms;
param.RSS.nMinFeatsPerc = nMinFeatsPerc;
param.RSS.nMaxFeatsPerc = nMaxFeatsPerc;

if act ~= 4 
    param = nk_RSS_config(param, setupfl, defaultsfl);        
end