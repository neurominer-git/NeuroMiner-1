% ********************************************************************
% * mex File compiling code for Random Forest (for linux)
% * mex interface to Andy Liaw et al.'s C code (used in R package randomForest)
% * Added by Abhishek Jaiantilal ( abhishek.jaiantilal@colorado.edu )
% * License: GPLv2
% * Version: 0.02
% ********************************************************************/
function compile_linux

    % if you have problems with the 
    fprintf('if you have problems with compilation using this method, use the alternative by uncommenting some lines in compile_linux.m');
    system('rm *.mexglx *.mexa64;');
    %system('make clean;make mex;');
    


    %Alternative compilation, comment out the system lines above and uncomment the below lines
    %that start with the mex and system command (total there are 3 lines to uncomment).
    system('gfortran -O2 -fpic -march=native -c src/rfsub.f -o rfsub.o')

%%sometimes you have to add -lgfortran at the end for both command lines below
    mex src/mex_ClassificationRF_train.cpp  src/classRF.cpp src/classTree.cpp src/rfutils.cpp src/cokus.cpp rfsub.o -output mexClassRF_train -lm -DMATLAB -O -v -lgfortran 
    mex src/mex_ClassificationRF_predict.cpp  src/classRF.cpp src/classTree.cpp src/rfutils.cpp src/cokus.cpp rfsub.o -output mexClassRF_predict -lm -DMATLAB -O -v