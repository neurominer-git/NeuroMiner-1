function nk_make(srcdir)
% =========================================================================
% function nk_make(srcdir)
% -------------------------------------------------------------------------
% Inputs:
% srcdir            :       source directory of files (<nm root>/cfiles)
% This function compiles all C files needed for NM
% =========================================================================
% (c) Nikolaos Koutsouleris 01/2017

cd(srcdir)

    make('sumskipnan.cpp')
    make('fastcorr.cpp')
    make('extrX.cpp')
    make('extrY.cpp')
    make('fdr.cpp')
    make('fdr_1d.cpp')
    make('ftest.cpp')
    make('intpow.cpp')
    make('mean_1d.cpp')
    make('mean_2d.cpp')
    make('pearson.cpp')
    make('spatconst27.cpp')
    make('symbolize.cpp')
    make('enumelem.cpp')
    mex ('fastAUC.cpp')

end

function make(funcname)

try
    fprintf('\nCompiling: %s\n',funcname)
    mex(funcname);
catch
    cprintf('*red','\nCompiling of %s failed!',funcname)
end

end
