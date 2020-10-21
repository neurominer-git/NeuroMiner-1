function rX = ranking(X)
% Transform an arbitrarily distributed signal to a uniformly distributed
% one.

nX = size(X,2);
sortX = zeros(size(X));
sI = sortX;
rX = sortX;

for n=1:nX
% Determine sort index 
    [sortX(:,n),sI(:,n)]=sort(X(:,n));
    % Resort ascending sequence according to 
    % the corresponding sort index. 
    rX(sI(:,n),n)= linspace(0.0,1.0,length(X(:,n)));
end 

return