function r = spear(x,y)

% Find the data length
N = length(x);

% Get the ranks of x
R = crank(x)';

% Get the ranks of y
S = crank(y)';
    
% Calculate the Spearman correlation coefficient
r = 1-6*sum((R-S).^2)/N/(N^2-1);

function r=crank(x)

u = unique(x); lU = length(u);
[~,z1] = sort(x);
[~,z2] = sort(z1);
r = (1:length(x))';
r = r(z2);

for i=1:lU;
    s = u(i) == x;
    r(s,1) = mean(r(s));
end






