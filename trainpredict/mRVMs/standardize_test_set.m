% authors: Ioannis Psorakis psorakis@gmail.com, Theodoros Damoulas
% theo@dcs.gla.ac.uk
%    Copyright NCR, 2009 
%
%    Use of this code is subject to the Terms and Conditions of the Research License Agreement,
%    agreed to when this code was downloaded, and a copy of which is available at 
%    http://www.dcs.gla.ac.uk/inference/pMKL/Terms_and_Conditions.html.

function Xtest_standardized = standardize_test_set(Xtest,Xtrain)

%sizeXtrain = size(Xtrain,1);
%sizeXtest = size(Xtest,1);

[Xtrain_standardized Xtrain_mean Xtrain_std] = standardize_data(Xtrain);

Xtest_standardized = standardize_data(Xtest,Xtrain_mean,Xtrain_std);

% XY = standardize_data([Xtrain;Xtest]);
% Ys = XY([end-sizeXtest+1:1:end],:);