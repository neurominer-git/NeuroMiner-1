function x_new = NextSolution2(x_curr, itt, n, featureID_order)
% Unlike NextSolution1, this routine NextSolution2 can pick which new voxel
% should or should not be included to the SA by considering the probability
% derived from MI-descending rank   
moduloItt = mod(itt,n); if moduloItt == 0, moduloItt = n; end
x_new = x_curr;
% @#$% user-defined
% ====================================================
% Approach II: introduce only some voxels
Delta = 2*((n-1)-(featureID_order(moduloItt)-1))/(n-1)-1;
x_new(featureID_order(moduloItt)) = rand(1) < 0.5+0.3*Delta; 
% ====================================================
x_new = logical(x_new);