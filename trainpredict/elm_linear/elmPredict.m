function scores = elmPredict( X, inW, bias, outW )
% FUNCTION predicts the labels of some testing data using a trained Extreme
% Learning Machine.
%
%	scores = elmPredict( X, inW, bias, outW );
%
% INPUT :
%	X				- data patterns (column vectors)
%	inW				- input weights vector (trained model)
%	bias			- bias vector (trained model)
%	outW			- output weights vector (trained model)
%
% OUTPUT :
%	scores			- prediction scores on the data
%

% number of test patterns
nTestData = size( X, 2 );

% compute the pre-H matrix
preH = inW * X;

% build the bias matrix
biasM = repmat( bias, 1, nTestData );

% update the pre-H matrix
preH = preH + biasM;

% apply the activation function
H = 1 ./ (1 + exp(-preH));

% compute prediction scores
scores = (H' * outW)';
