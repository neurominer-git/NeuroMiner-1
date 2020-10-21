function matLearn = nk_matLearn_config(matLearn, framework, act)

if ~exist('framework','var') || isempty(framework)
    matLearn.learner.framework = char(nk_input('Select the learning framework',0,'mq', ...
        'Regression|Binary classification|Multi-group classification|Ordinal regression|Multi-label prediction|Unsupervised learning', ...
        {'regression','binaryclass','multiclass','ordinal','multilabel','unsupervised'},1));
else
     matLearn.learner.framework = framework;
end
if ~isfield(matLearn,'learner'), matLearn.learner = []; end
if ~isfield(matLearn.learner,'options'), matLearn.learner.options = []; end
matLearn.learner = nk_matLearn_getopts_config(matLearn.learner,'get_learners',[],matLearn, matLearn.learner.framework);

nk_PrintLogo
if ~exist('act','var') || isempty(act)
    act = nk_input('Select action',0,'mq','Select matLearn algorithm|Define respective parameters',[1,2],1);
    flg = true;
else
    flg = false;
end

switch act
    case 0
        return
    case 1
        matLearn = select_algo(matLearn);
    case 2
        matLearn.learner.options = nk_matLearn_getopts_config(matLearn.learner.options,'get_learner_params', char(matLearn.algo), matLearn, matLearn.learner.framework);
        if ~isempty(matLearn.learner.options)
            matLearn = nk_matLearn_DefAlgoParams_config(matLearn, matLearn.learner.options, true);
        else
            matLearn.Params = [];
        end
    case 3
        matLearn.learner.options = nk_matLearn_getopts_config(matLearn.learner.options,'get_learner_params', char(matLearn.algo), [], matLearn.learner.framework);
        if ~isempty(matLearn.learner.options)
            matLearn = nk_matLearn_DefAlgoParams_config(matLearn, matLearn.learner.options, false);
        else
            matLearn.Params = [];
        end
end

if flg 
    matLearn = nk_matLearn_config(matLearn, matLearn.learner.framework);
end

% ==========================================================================
function param = select_algo(param)

param.algo = nk_matLearn_IO_config(param.learner, param, 1);

