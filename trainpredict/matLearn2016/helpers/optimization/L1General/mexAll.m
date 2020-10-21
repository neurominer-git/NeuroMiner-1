% minFunc
fprintf('Compiling minFunc files...\n');
mex minFunc/mcholC.c
mex minFunc/lbfgsC.c