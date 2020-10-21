function [inW, bias, outW, scores] = elmTrain( X, Y, nHiddenNeurons, C )
% FUNCTION trains the Extreme Learning Machine. The activation function is
% sigmoid which could be changed easily if needed.
%
%	[inW, bias, outW, scores] = elmTrain( X, Y, nHiddenNeurons, C );
%
% INPUT :
%	X				- data patterns (column vectors)
%	Y				- numeric labels for each pattern (1, ... )
%	nHiddenNeurons	- number of hidden neurons
%	C				- regularization parameter
%
% OUTPUT :
%	inW				- input weights matrix
%	bias			- bias vector
%	outW			- output weights matrix
%	scores			- scores on the own training data
%

nTrainData = size( X, 2 );
nInputNeurons = size( X, 1 );
nClasses = length( unique(Y) );

targets = zeros( nClasses, nTrainData );

% populate 1-of-k target matrix - transforming labels vector
for i = 1 : nTrainData
	targets( Y(i), i ) = 1;
end
targets = targets * 2 - 1;	% from 0/1 to -1/1

% generate random input weight matrix
inW = rand( nHiddenNeurons, nInputNeurons ) * 2 - 1;

% generate random hidden neuron vector
bias = rand( nHiddenNeurons, 1 );

% compute the pre-H matrix
preH = inW * X;

% build the bias matrix
biasM = repmat( bias, 1, nTrainData );

% update the H matrix
preH = preH + biasM;

% calculate hidden neuron output matrix H
H = 1 ./ (1 + exp(-preH));

% compute regularized output weights matrix
outW = ( eye(nHiddenNeurons)/C + H * H') \ H * targets';

% output for the training data
scores = (H' * outW)';
