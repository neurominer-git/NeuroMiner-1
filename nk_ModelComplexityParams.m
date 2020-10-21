function [sumNSV, meanNSV, stdNSV, rSV] = nk_ModelComplexityParams(M, S)
% -------------------------------------------------------------------------
% function [sumNSV, meanNSV, stdNSV] = nk_ModelComplexityParams(M)
% -------------------------------------------------------------------------
% Compute model complexity depending on the underlying model type.
% 
% Input:
% M :       Models to be analyzed
% S :       Number of cases in fold CV1 cycle
%
% Output:
% sumNSV :  sum of support vectors (SV)      
% meanNSV : mean number of SV       
% stdNSV :  std number of SV
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris 07/2020

global SVM

[ix,jx,nclass] = size(M); 
if strcmp(SVM.prog,'SEQOPT') 
    if ~isfield(SVM.SEQOPT,'complexity_mode') || isempty(SVM.SEQOPT.complexity_mode) 
        COMPLEXITY = 1;
    else
        COMPLEXITY = SVM.SEQOPT.complexity_mode;
    end
end

sumNSV = zeros(nclass,1);
meanNSV = zeros(nclass,1);
stdNSV = zeros(nclass,1);
rSV = zeros(nclass,1);

for curclass=1:nclass
    
    NSV=[];
    switch SVM.prog

        case {'LIBSVM','MikRVM','MVTRVR','MSTOOL'}
            for i=1:ix
                for j=1:jx
                    for l=1:length(M{i,j,curclass})
                        NSV = [NSV; M{i,j,curclass}{l}.totalSV];
                    end
                end
            end
        case 'MKLRVM'
            for i=1:ix
                for j=1:jx
                    for l=1:length(M{i,j,curclass})
                        NSV = [NSV; M{i,j,curclass}{l}.N_prototypical];
                    end
                end
            end
        case 'SEQOPT'
            for i=1:ix
                for j=1:jx
                    for l=1:length(M{i,j,curclass})
                        % Currently we can compute 'complexity' of the 
                        % sequential prognostication algorithm ... 
                        switch COMPLEXITY
                            case 1
                                % ... as the percentage of cases that reach 
                                % the last examination node.
                                mX = M{i,j,curclass}{l}.examsfreq(end);
                            case 2
                                % ... as the average percentage of cases
                                % propagated across the chain / number of
                                % nodes in the chain. Problem of this
                                % approach is that sequences with only one
                                % node will be heavily favored if
                                % complexity is used to regularize
                                % performance
                                mX = sum( 100 - M{i,j,curclass}{l}.examsfreq(1:end-1) ) ...
                                    / size(SVM.SEQOPT.C,1);
                        end
                        NSV = [NSV; mX];
                    end
                end
            end
        otherwise
            NSV = S(curclass);
    end

    sumNSV(curclass)  = sum(NSV);
    meanNSV(curclass) = mean(NSV);
    stdNSV(curclass)  = std(NSV);
    switch SVM.prog
        case 'SEQOPT'
            rSV(curclass)     = meanNSV(curclass);
        otherwise
            rSV(curclass)     = meanNSV(curclass)/S(curclass)*100;
    end
end

