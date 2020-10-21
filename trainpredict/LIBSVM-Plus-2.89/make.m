% This make.m is used under Windows

% add -largeArrayDims on 64-bit machines

mex -largeArrayDims -O -c svm.cpp
mex -largeArrayDims -O -c svm_model_matlab.c
mex -largeArrayDims -O svmtrain289PLUS.c svm.obj svm_model_matlab.obj
mex -largeArrayDims -O svmpredict289PLUS.c svm.obj svm_model_matlab.obj
%mex -O libsvmread.c
%mex -O libsvmwrite.c
