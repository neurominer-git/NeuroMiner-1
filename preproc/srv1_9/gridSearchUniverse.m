function [alphaOptimal,betaOptimal,accMax]=gridSearchUniverse(trainSet,trainClass,option)
% Line or Grid Search
%%%%
% Copyright (C) <2012>  <Yifeng Li>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Contact Information:
% Yifeng Li
% University of Windsor
% li11112c@uwindsor.ca; yifeng.li.cn@gmail.com
% July 04, 2012
%%%%

optionDefault.classifier='svmLi';
optionDefault.alphaRange=2:5; % exponent of 2
optionDefault.betaRange=6:10; % exponent of 2
optionDefault.kfold=3;
optionDefault.rerun=1;
%optionDefault.normalization=0;
optionDefault.optCl=[];

if nargin==2
    option=optionDefault;
else
    option=mergeOption(option,optionDefault);
end
numClass=numel(unique(trainClass));
numAlphaRange=numel(option.alphaRange);
numBetaRange=numel(option.betaRange);
AllAccs=zeros(numAlphaRange,numBetaRange);

for c=1:numAlphaRange
    for i=1:numBetaRange
        switch option.classifier
            case {'svm','svmLi','lsvmLi','svmm'}
                alpha=(option.alphaRange(c));
                beta=2^(option.betaRange(i));
                option.optCl.C=alpha;
                option.optCl.ctrl=alpha;
                option.optCl.param=beta;
            case 'gdtree'
                alpha=2^(option.alphaRange(c));
                beta=2^(option.betaRange(i));
                option.optCl.optCl.C=alpha;
                option.optCl.optCl.ctrl=alpha;
                option.optCl.optCl.param=beta;
            case {'nnls','src2'}
                alpha=2^(option.alphaRange(c));
                beta=2^(option.betaRange(i));
                option.optCl.lambda=alpha;
                option.optCl.param=beta;
            case 'hdlm'
                alpha=2^(option.alphaRange(c));
                beta=2^(option.betaRange(i));
                if ~strcmp(option.optCl.kernel,'sigmoid')
                    alpha=[];
                end
                option.optCl.param=[alpha,beta];
            case 'svdd'
                alpha=option.alphaRange(c);
                beta=2^(option.betaRange(i));
                option.optCl.nu=alpha;
                option.optCl.param=beta;
            case 'svdda'
                alpha=option.alphaRange(c);
                beta=2^(option.betaRange(i));
                option.optCl.fracsv=alpha;
                option.optCl.param=beta;
            case 'svddm'
                alpha=option.alphaRange(c);
                beta=option.betaRange(i);
                option.optCl.nu=alpha;
                option.optCl.nus=beta;
%                 option.optCl.param=2^1;
            otherwise
                error('Please select a correct classifier!');
        end
        perf=cvExperiment(trainSet,trainClass,option.kfold,option.rerun,'none',[],'none',[],{option.classifier},{option.optCl});
        AllAccs(c,i)=perf(end-1);
    end
end
[AllAccsRow,indsRow]=max(AllAccs,[],1);
[accMax,indCol]=max(AllAccsRow);
indRow=indsRow(indCol);
betaOptimal=option.betaRange(indCol);
alphaOptimal=option.alphaRange(indRow);
end