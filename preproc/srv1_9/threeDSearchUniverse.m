function [alphaOptimal,betaOptimal,gammaOptimal,bestAcc,AllAccs]=threeDSearchUniverse(trainSet,trainClass,option)
% Three Dimension Search
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
% September 12, 2012
%%%%

optionDefault.classifier='svddm';
optionDefault.alphaRange=2:5; % may be exponent of 2
optionDefault.betaRange=2:5; % may be exponent of 2
optionDefault.gammaRange=2:5; % may be exponent of 2
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
numGammaRange=numel(option.gammaRange);
AllAccs=zeros(numAlphaRange,numBetaRange,numGammaRange);
bestAcc=-1;
alphaOptimal=NaN;
betaOptimal=NaN;
gammaOptimal=NaN;
for c=1:numAlphaRange
    for i=1:numBetaRange
        for k=1:numGammaRange
            switch option.classifier
                case {'svm','svmLi','lsvmLi','svmm'}
                    alpha=(option.alphaRange(c));% 2^(option.alphaRange(c));
                    beta=2^(option.betaRange(i));
                    gamma=2^(option.gammaRange(k));
                    if ~strcmp(option.optCl.kernel,'sigmoid')
                        gamma=[];
                    end
                    option.optCl.C=alpha;
                    option.optCl.ctrl=alpha;
                    option.optCl.param=[beta;gamma];
                case 'gdtree'
                    alpha=2^(option.alphaRange(c));
                    beta=2^(option.betaRange(i));
                    gamma=2^(option.gammaRange(k));
                    if ~strcmp(option.optCl.kernel,'sigmoid')
                        gamma=[];
                    end
                    option.optCl.optCl.C=alpha;
                    option.optCl.optCl.ctrl=alpha;
                    option.optCl.optCl.param=[beta;gamma];
                case {'hdlm','nnls','nc','lrc'}
                    alpha=2^(option.alphaRange(c));
                    beta=2^(option.betaRange(i));
                    gamma=2^(option.gammaRange(k));
                    option.optCl.param=[alpha;beta;gamma];
                case 'svdd'
                    alpha=option.alphaRange(c);
                    beta=2^(option.betaRange(i)); 
                    gamma=2^(option.gammaRange(k));
                    option.optCl.nu=alpha;
                    option.optCl.param=[beta;gamma];
                case 'svddm'
                    alpha=option.alphaRange(c);
                    beta=option.betaRange(i);
                    gamma=2^(option.gammaRange(k));
                    option.optCl.nu=alpha;
                    option.optCl.nus=beta;
                    option.optCl.param=gamma;
                otherwise
                    error('Please select a correct classifier!');
            end
            perf=cvExperiment(trainSet,trainClass,option.kfold,option.rerun,'none',[],'none',[],{option.classifier},{option.optCl});
            if perf(end-1)>bestAcc
                bestAcc=perf(end-1);
                alphaOptimal=alpha;
                betaOptimal=beta;
                gammaOptimal=gamma;
            end
            AllAccs(c,i,k)=perf(end-1);
        end
    end
end
end