function [w,L] = sbmlr(x, t, w)
%
% SBMLR - sparse Bayesian multinomial logistic regression 
%
%    ALPHA = SBMLR(X,T) computes the weights for a sparse multinomial
%    regression model, where the sparsity is obtained using Bayesian
%    regularisation with a Laplace prior, as described in [1].  Here X is
%    a matrix of input features, such that each column represents a feature
%    and each row represents a pattern.  T is a matrix of target patterns,
%    such that each column represents a class and each row represents a
%    pattern.  All of the elements of T are assumed to lie in the interval
%    (0,1) and that the rows of T sum to one.
%
%    References
%
%    [1] G. C. Cawley, N. L. C. Talbot and M. Girolami, "Sparse multinomial
%        logistic regression via Bayesian regularisation using a Laplace
%        prior", Neural Information Processing Systems 19, 2006.
%
%    [4] P. M. Williams, "Bayesian regularisation and pruning using a Laplace
%        prior", Neural Computation, 7(1):117-143, 1995.
%
%    [5] W. L. Buntine and A. S. Weigend, "Bayesian back-propagation", Complex
%        Systems, 5:603-643, 1991.

%
% File        : sbmlr.m
%
% Date        : Monday 17th April 2006
%
% Author      : Gavin C. Cawley
%
% Description : Sparse multinomial logistic regression using a Laplace prior.
%
% History     : 17/04/2006 - v1.00
%
% Copyright   : (c) Dr Gavin C. Cawley, April 2006.
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
%

if nargin == 0

   w.MaxEpoch  = 100000;
   w.MaxIter   = 100;
   w.Epsilon   = 1e-3;
   w.Pi        = 1e-1;
   w.Lambda    = 1e+0;
   w.LambdaMax = realmax;
   w.LambdaMin = realmin;
   w.Tau       = 1e-9;
   w.TolFun    = 1e-6;

   return;

else

   opts      = sbmlr;
   [ntp,d]   = size(x);
   c         = size(t,2);
   preselect = [];

   if nargin == 2

      w     = zeros(d, c);
      Ew    = 0.0;
      alpha = 0.0;

   else

      Ew    = sum(abs(w(:)));

      if Ew > eps

         alpha = sum(w(:)~=0)/Ew;

      else

         alpha = 0.0;

      end

   end

end

deleted = zeros(size(w));
turnip  = 1.0;
exp_a   = exp(x*w);
y       = exp_a./repmat(sum(exp_a,2), 1, size(w,2));
Ed      = -sum(log(y(t(:)==1) + realmin));

for i=1:opts.MaxEpoch

   W           = sum(w(:) ~= 0);
   E_w         = sum(abs(w(w~=0)));
   L           = Ed + alpha*Ew;
   Y           = min(max(y, 1-1e6), 1-1e-6);
   dEd         = ((Y - t)'*x)';
   d2Ed        = (Y'.*(1 - Y)'*(x.^2))';
   grad        = zeros(size(w));
   idx         = find(w > 0.0);
   grad(idx)   = dEd(idx) + alpha;
   idx         = find(w < 0.0);
   grad(idx)   = dEd(idx) - alpha;
   idx         = find((w == 0.0) & ((dEd + alpha) < 0.0));
   grad(idx)   = dEd(idx) + alpha;
   idx         = find((w == 0.0) & ((dEd - alpha) > 0.0));
   grad(idx)   = dEd(idx) - alpha;
   delta       = abs(grad.*grad)./d2Ed; %./d2Ed; %1.5*grad.*grad./d2Ed;
   idx         = find(deleted > 10);
   delta(idx)  = 0;
   [delta,idx] = max(delta(:)); 

   fprintf(1, 'epoch = %-4d : alpha = %f : L = %f : W = %d     \r', i, alpha, L, W);

   if delta > 0

      [feature,class] = ind2sub(size(w), idx);

      if w(idx) > 0.0

         weight = w(idx);
         scale  = 1.0;

         for k=1:10
   
            w(idx)         = max(weight - scale*grad(idx)/d2Ed(idx), 0.0);
            W              = sum(w(:) ~= 0);
            exp_a(:,class) = exp(x*w(:,class));
            y              = exp_a./repmat(sum(exp_a,2), 1, size(w,2));
            Ed             = -sum(log(y(t(:)==1) + realmin));
            Ew             = sum(abs(w(:)));
            L_new          = Ed + alpha*Ew;

            if L_new > L

               scale = 0.5*scale;

            else

               break;

            end

         end

         if w(idx) == 0.0

            deleted(idx) = deleted(idx) + 1;

         end

      elseif w(idx) < 0.0

         weight = w(idx);
         scale  = 1.0;

         for k=1:10

            w(idx)         = min(weight - scale*grad(idx)/d2Ed(idx), 0.0);
            W              = sum(w(:) ~= 0);
            exp_a(:,class) = exp(x*w(:,class));
            y              = exp_a./repmat(sum(exp_a,2), 1, size(w,2));
            Ed             = -sum(log(y(t(:)==1) + realmin));
            Ew             = sum(abs(w(:)));
            L_new          = Ed + alpha*Ew;

            if L_new > L

               scale = 0.5*scale;

            else

               break;

            end

         end

         if w(idx) == 0.0

            deleted(idx) = deleted(idx) + 1;

         end

      else

         weight = w(idx);
         scale  = 1.0;

         for k=1:10

            w(idx)         = weight - scale*grad(idx)/d2Ed(idx);
            W              = sum(w(:) ~= 0);
            exp_a(:,class) = exp(x*w(:,class));
            y              = exp_a./repmat(sum(exp_a,2), 1, size(w,2));
            Ed             = -sum(log(y(t(:)==1) + realmin));
            Ew             = sum(abs(w(:)));
            L_new          = Ed + alpha*Ew;

            if L_new > L

               scale = 0.5*scale;

            else

               break;

            end

         end

      end

      if Ew > 1e-6

         alpha = sum(w(:)~=0)/Ew;

      end

      turnip = scale;

      if abs((L_new - L)/L_new) < opts.TolFun

         break;

      end

   else

      % can no longer make progress

      break;

   end

end

% compute marginal likelihood

idx             = find(w(:)~=0);
[feature,class] = ind2sub(size(w), idx);

for i=1:length(idx)

   for j=1:length(idx)

      if class(i) == class(j)

         A(i,j) = ( y(:,class(i)) - y(:,class(i)).*y(:,class(j)) )'*( x(:,feature(i)).*x(:,feature(j)) );

      else

         A(i,j) = -(y(:,class(i)).*y(:,class(j)))'*(x(:,feature(i)).*x(:,feature(j)));

      end

   end

end

if W == 0

   L = Ed;

else

   L = Ed + W*log(sum(abs(w(:)))) - gammaln(W) + W*log(2) + 0.5*log(det(A));

end

fprintf(1, '                                                                            \r');

% bye bye...

