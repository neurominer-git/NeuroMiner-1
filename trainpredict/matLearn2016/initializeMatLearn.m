function initializeMatLearn(matLearnDir)
%% Run this script to initialize matLearn

disp('Initializing MatLearn...');

%% change Matlab's working directory to this script's directory
addpath(matLearnDir);

%% include all MatLearn files in path
addpath(genpath(matLearnDir));

end
