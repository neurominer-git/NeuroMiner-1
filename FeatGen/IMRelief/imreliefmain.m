% =========================================================================
% FORMAT [R, nsigma, nlambda] = imreliefmain(Y, labels, SortInd)
% =========================================================================
% This function provides an interface between NeuroMiner and the IMRelief
% feature selection function. It invokes either the CPU or the GPU version
% of IMRelief depending on the hardware / user settings and performs
% multiple invocations of IMRelief if either lambda or sigma parameters
% ranges have been defined by the used
%
% Inputs:
% Y :           [m x n] matrix, with m samples / patterns and n features
% labels :      m x 1 vector containing supervision info (currently only
%               two-class problems supported)
% (SortInd :    m x 1 resampling vector)
%
% Outputs:      
% R :           n x 1 Weight vector
% nsigma :      number of sigma parameters tested
% nlambda :     number of lambda parameters tested
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 07/2011
function [R, nsigma, nlambda] = imreliefmain(Y, labels, curclass, SortInd)

global FEATSEL CUDA CPUAVAIL BATCH

if exist('SortInd','var') && ~isempty(SortInd), Y = resamp(Y, SortInd, CPUAVAIL); end;

if BATCH, verbose = false; else verbose = true; end

if iscell(FEATSEL.imrelief)
    if numel(FEATSEL.imrelief)<curclass
        imrelief = FEATSEL.imrelief{end};
    else
        imrelief = FEATSEL.imrelief{curclass};
    end
else
    imrelief = FEATSEL.imrelief;
end

sigma = imrelief.sigma;
lambda = imrelief.lambda;
distance = imrelief.distance;
nsigma = numel(imrelief.sigma);
nlambda = numel(imrelief.lambda);
plotfigure = imrelief.plotfigure;
maxiter = imrelief.maxiter;

if nsigma > 1 || nlambda > 1, end
if size(Y,2) ~= length(labels), Y = Y'; end
if ~CUDA
    IMRelief = 'IMRelief_Sigmoid_FastImple';
else
    IMRelief = 'IMRelief_Sigmoid_FastImple_gpu';
    Y = GPUsingle(Y);
end

% Eliminate useless features
[ Y, NonPruneVec ] = nk_PerfElimZeroObj(Y', []); Y=Y';

if CPUAVAIL > 1 && ~CUDA && nsigma > 1
    if CPUAVAIL <= nsigma,
        cpus = floor(nsigma/CPUAVAIL);
    else
        cpus = nsigma;
    end
    R = cell(nsigma, nlambda);
    feval('matlabpool','open',cpus)
    parfor i=1:nsigma
        for j=1:nlambda
            R{i}{j} = zeros(numel(NonPruneVec),1);
            DScore = feval(IMRelief, Y, labels, ...
                                    distance, ...
                                    sigma(i), ...
                                    lambda(j), ...
                                    maxiter, ...
                                    plotfigure, verbose);
            R{i}{j}(NonPruneVec) = DScore;
        end
    end
    matlabpool close 
    R = cell2mat(R);
    size(R)
else
    R = zeros(numel(NonPruneVec),nsigma, nlambda);
    for i=1:nsigma
        for j=1:nlambda
            DScore = feval(IMRelief, Y, labels, ...
                                    distance, ...
                                    sigma(i), ...
                                    lambda(j), ...
                                    maxiter, ...
                                    plotfigure, verbose);
            if CUDA, 
                R(NonPruneVec,i,j) = double(DScore);
            else
                R(NonPruneVec,i,j) = DScore;
            end
        end
    end

end
