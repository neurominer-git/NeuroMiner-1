function [sY, IN] = nk_PerfDiscretizeObj(Y, IN)
% =========================================================================
% FORMAT [dY, IN] = nk_PerfDiscretizeObj(Y, IN)
% =========================================================================
% Do one of the following to digitize feature matrix Y
%
% Either:
% Discretize Y columnwise to mean +/- alpha*std, where alpha is a
% vector of discrete values. By default alpha is set to 0:0.5:4.
% 
% Or:
% Symbolize Y columnwise using an entropy-based encoding of features into
% discrete code series. Thus, features with high entropy may have a higher
% number of discrete values than features with low entropy.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris 02/2017    

% =========================== WRAPPER FUNCTION ============================ 

if iscell(Y) && exist('IN','var') && ~isempty(IN)
    sY = cell(1,numel(Y)); 
    for i=1:numel(Y), [sY{i}, IN] =  IN.method(Y{i}, IN); end
else
    [ sY, IN ] = IN.method(Y, IN );
end
% =========================================================================
