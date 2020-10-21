function [X,Y] = proj_def(X, Y, u, v)
%
%   Perform projection deflation of matrices X and Y using u and v
%   Please check Monteiro et al. 2016 for details:
%   doi:10.1016/j.jneumeth.2016.06.011
%
%   Inputs: X, Y    - data matrices in the form: samples x features.
%
%           u, v    - weight vectors for X and Y, respectively
%
%
%   Outputs: X, Y    - deflated data matrices
%
%
%   Version: 2016-07-06
%__________________________________________________________________________

% Written by Joao Matos Monteiro
% Email: joao.monteiro@ucl.ac.uk

X = X - (X*u)*u';
Y = Y - (Y*v)*v';

end