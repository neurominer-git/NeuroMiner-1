function x_new = NextSolution1(x_curr, itt, n, featureID_order)
% The new voxel is introduced to the solution in the order determined by
% futureID_order.
moduloItt = mod(itt,n); if moduloItt == 0, moduloItt = n; end
x_new = x_curr;
% @#$% user-defined
% ====================================================
% Appproach I: introduce every voxel
x_new(featureID_order(moduloItt)) = 1; 
% ====================================================
x_new = logical(x_new);