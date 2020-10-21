/******************************************************************************

File        : blogreg.c

Date        : Thursday 30th March 2006

Author      : Dr Gavin C. Cawley

Description : MEX implementation of a Bayesian logsitic regression method
              based on the training algorithm for Shevade and Keerthi's sparse
	      logistic regression procedure [1,2] (see blogreg.m for usage).
              Details of the algorithm are given in [3].

References  : [1] Shevade, S. K. and Keerthi, S. S., "A simple and effecient 
                  algorithm for gene selection using sparse logistic
		  regression", Technical Report CD-02-22, Control Division,
		  Department of Mechanical Engineering, National University
		  of Singapore, Singapore - 117 576, 2002.

              [2] Shevade, S. K. and Keerthi, S. S., "A simple and effecient 
                  algorithm for gene selection using sparse logistic
		  regression", Bioinformatics, vol. 19, no. 17, pp 2246-2253,
		  2003.

              [3] Cawley, G. C. and Talbot, N. L. C., "Gene selection in
                  cancer classification using sparse logistic regression with
                  Bayesian regularisation", Bioinformatics (submitted), 2006.

History     : 02/10/2004 - v1.00
              30/03/2006 - v1.10 minor improvements to comments etc.

Copyright   : (c) G. C. Cawley, March 2006

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

******************************************************************************/

#include <float.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>

#include "mex.h"

/* symboloc constants */

const int I_z  = 0;
const int I_nz = 1;

/* global data */

int     *set;         
double  *alpha;
double  *xi;
double  *exp_xi;
double   lambda;      /* regularisation parameter        */
int      ntp;         /* number of training patterns     */
int      d;           /* number of input features        */
double  *F;
double  *delta;
double  *tmp1;
double  *tmp2;
double **x;
double  *t;
double   tol;
double   E_alpha;     /* sum of magnitude of weights     */
double   E_d;         /* the data misfit term            */
int      N;           /* number of active features       */
double   M;           /* the standard objective function */
double   Q;           /* the Bayesian objective function */
int      epoch;

/* macros */

#define SWAP(a,b,type) { auto type a_b_type = a ; a = b ; b = a_b_type; }

/* utility functions */

double square(double x)
{
   return x*x;
}

double doubleMax3(double a, double b, double c)
{
   return a > b ? (a > c ? a : c) : (b > c ? b : c);
}

/******************************************************************************

Procedure   : objective 

Parameters  : none

Returns     : double - value of objective function

Description : Evaluate the objective function for the current set of
              coefficients.

******************************************************************************/

void compute_objective()
{
   int i;

   /* evaluate the data misfit term */

   E_d = 0.0;

   for (i = 0; i < ntp; i++)
   {
      E_d += log(1.0 + exp_xi[i]); 
   }

   /* evaluate objective functions */

   M = E_d + lambda*E_alpha;
   Q = E_d + (E_alpha == 0.0 ? 0.0 : N*log(E_alpha));
}

/******************************************************************************

Procedure   : updateFsAndDeltas

Parameters  : int s - the set of input features to update (I_z or I_nz)

Returns     : void

Description : Update the F and delta statistics for all input features
              beling to set I_z or I_nz.

******************************************************************************/

void updateFsAndDeltas(int s)
{
   int i, j;

   double delt, f;

   for (i = 0; i < d; i++)
   {
      if (set[i] == s) 
      {
         f    = 0.0;
         delt = 0.0;

         for (j = 0; j < ntp; j++)
         {
            f    += exp_xi[j]*t[j]*x[i][j]/(1 + exp_xi[j]);
            delt += exp_xi[j]*square(x[i][j])/square(1 + exp_xi[j]);
         }

         F[i]     = f;
         delta[i] = delt;
      }
   }
}

/******************************************************************************

Procedure   : updateFAndDelta

Parameters  : int feature

Returns     : void

Description : Update the F and delta statistics for a specified input feature.

******************************************************************************/

void updateFAndDelta(int feature)
{
   int i;

   double delt = 0.0, f = 0.0;

   for (i = 0; i < ntp; i++)
   {
      f    += exp_xi[i]*t[i]*x[feature][i]/(1 + exp_xi[i]);
      delt += exp_xi[i]*square(x[feature][i])/square(1 + exp_xi[i]);
   }

   F[feature]     = f;
   delta[feature] = delt;
}

/******************************************************************************

Procedure   : updateXi

Parameters  : int    feature - the feature with the most recently updated
                               coefficient
              double change  - size of the change in the coefficient

Returns     : void

Description : Update the xi and exp_xi statistics for all features, following
              a given change in the coefficient for specified feature.

******************************************************************************/

void updateXi(int feature, double change)
{
   int i;

   for (i = 0; i < ntp; i++)
   {
      xi[i] += t[i]*change*x[feature][i];

      exp_xi[i] = exp(xi[i]);
   }
}

/******************************************************************************

Procedure   : findMaxViolator

Parameters  : int s - the set of input features to search (I_z or I_nz)

Returns     : int - index of the maximally violating coefficient

Description : Return the index of the coefficient that maximally violates
              the optimality conditions from the specified set.  Returns -1
	      if no feature violates the optimality conditions by more than
	      the tolerance.

******************************************************************************/

int findMaxViolator(int s)
{
   int i, violator;

   double maxViol, viol;

   if (s == I_nz)
   {
      updateFsAndDeltas(I_nz);

      if (fabs(F[0]) > tol)
      {
         maxViol  = fabs(F[0]);
         violator = 0;
      }
      else
      {
         maxViol  = tol;
	 violator = -1;
      }


      for (i = 1; i < d; i++)
      {
         if (set[i] == I_nz)
	 {
            viol = fabs(alpha[i] > 0.0 ? lambda - F[i] : lambda + F[i]);

	    if (viol > maxViol)
	    {
               maxViol  = viol;
	       violator = i;
	    }
	 }
      }
   }
   else
   {
      double scale = N > 0 ? (N+1.0)/N : 1.0;

      updateFsAndDeltas(I_z);

      violator = -1;
      maxViol  = tol;

      for (i = 1; i < d; i++)
      {
         if (set[i] == I_z)
	 {
            viol = doubleMax3(F[i] - scale*lambda, -F[i] - scale*lambda, 0.0);

	    if (viol > maxViol)
	    {
               maxViol  = viol;
	       violator = i;
            }
	 }
      }
   }

   return violator;
}

/******************************************************************************

Procedure   : computeSlopeAtZero

Parameters  :

Returns     :

Description :

******************************************************************************/

double evaluateSlope(int feature)
{
   if (feature == 0)
   {
      return F[0];
   }
   else if (alpha[feature] < 0.0)
   {
      return -lambda - F[feature];
   }
   else if (alpha[feature] > 0.0)
   {
      return lambda - F[feature];
   }
   else if (alpha[feature] == 0.0)
   {
      if ((lambda - F[feature]) < 0.0)
      {
         return lambda - F[feature];
      }
      else if ((-lambda - F[feature]) > 0.0)
      {
         return -lambda - F[feature];
      }
      else
      {
         return 0.0;
      }
   }
}

/******************************************************************************

Procedure   : optimiseAlpha

Parameters  : int violator

Returns     : none

Description :

******************************************************************************/

void optimiseAlpha(int violator)
{
   int i, flag;

   double L = 0.0, H = 0.0, a, F_tmp, delta_tmp, slope, d_right, d_left;

   double aleph = alpha[violator];

   /* bracket minimum with interval excluding zero (except at a boundary) */

   if (violator == 0)
   {
      /* case 1 - unregularised bias term */

      L     = -DBL_MAX;
      H     = +DBL_MAX;
      slope = -F[0];
   }
   else
   {
      /* the violator is a regularised parameter */

      if (aleph == 0.0)
      {
         d_right = +lambda - F[violator];
         d_left  = -lambda - F[violator]; 

         if (d_right < 0.0)
         {
            /* case 2 - alpha = 0 and right derivative is negative */
   
            L     = 0.0;
	    H     = DBL_MAX;
	    slope = d_right; 
         }
         else if (d_left > 0.0)
         {
            /* case 3 - alpha = 0 and left derivative is positive */

            L     = -DBL_MAX;
	    H     = +0.0;
	    slope = d_left; 
         }
         else
         {
            /* case 4 - parameter stable at zero (should never happen) */

            set[violator] = I_z;
      
	    return;
         }
      }
      else if (aleph != 0.0)
      {
         slope = (aleph > 0.0 ? +lambda : -lambda) - F[violator];

         if (aleph > 0.0 && slope < 0.0)
         {
            /* case 5 - aleph positive and derivative is negative */
   
            L = aleph;
            H = +DBL_MAX;
         }
         else if (aleph < 0.0 && slope > 0)
         {
            /* case 6 - alpha negative and derivative is positve */

            L = -DBL_MAX;
            H = aleph;
         }
         else
         {
            /* store a copy of xi, exp_xi, F[violator] and delta[violator] */

            memcpy(tmp1,     xi, ntp*sizeof(double));
            memcpy(tmp2, exp_xi, ntp*sizeof(double));

            F_tmp     = F[violator];
            delta_tmp = delta[violator];

            /* update statistics assuming alpha[violator] = 0.0 */
   
            updateXi(violator, aleph);
            updateFAndDelta(violator);

	    alpha[violator] = 0.0;
	    set[violator]   = I_z;
	    N               = N - 1;
	    E_alpha         = E_alpha - fabs(aleph);

	    /* compute left and right derivaltives */

            d_right = +lambda - F[violator];
            d_left  = -lambda - F[violator]; 

	    if ((d_right>0.0 && d_left<0.0) || d_right==0.0 || d_left==0.0)
	    {
               /* parameter can be safely pruned */

	       return;
	    }
	    else
            {
               /* bracket minimum */

               if (aleph > 0.0 && slope > 0.0 && d_right > 0.0)
               {
                  /* case 2 */
   
	          aleph = 0.0;
	          slope = d_left;
	          L     = -DBL_MAX;
	          H     = +0.0;
               }
               else if (aleph < 0.0 && slope < 0.0 && d_left < 0.0)
               {
                  /* case 5 */
   
	          aleph = 0.0;
	          slope = d_right;
	          L     = 0.0;
	          H     = DBL_MAX;
               }
               else
	       {
	          if (aleph > 0.0 && slope > 0.0 && d_right < 0.0)
                  {
                     /* case 3 */
      
                     L = 0.0;
	             H = aleph;
	          }
                  else if (aleph < 0 && slope < 0 && d_left > 0)
                  {
                     /* case 6 */
      
	             L = aleph; 
	             H = 0.0;
                  }

                  /* restore xi and exp_xi etc. */
   
	          SWAP(    xi, tmp1, double*)
	          SWAP(exp_xi, tmp2, double*)
   
	          alpha[violator] = aleph;
	          set[violator]   = I_nz;
	          F[violator]     = F_tmp;
	          delta[violator] = delta_tmp;
	          N               = N + 1;
	          E_alpha         = E_alpha + fabs(aleph);
               }
            }
         }
      }
   }

   /* optimise coefficient for violator via Newton's method */

   for (i = 0; i < 100; i++)
   {
      if (fabs(slope) < 0.1*tol)
      {
         break;
      }

      if (i == 0)
      {
         a = aleph - 0.5*slope/(delta[violator]);
      }
      else
      {
         a = aleph - slope/(delta[violator]);
      }

      if (a < L || a > H)
      {
         a = 0.5*(L+H);
      }

      updateXi(violator, aleph - a);

      updateFAndDelta(violator);

      if (violator == 0)
      {
         slope = -F[violator];
      }
      else
      {
         if (a > 0.0)
	 {
            slope = lambda - F[violator];
	 }
	 else
	 {
            slope = -lambda - F[violator];
	 }
      }

      if (slope > 0.1*tol)
      {
         H = a;
      }
      else if (slope < -0.1*tol)
      {
         L = a;
      }

      aleph = a;
   }

   /* update coefficient vector */

   if (violator != 0)
   {
      N       = alpha[violator] == 0.0 && aleph != 0.0 ? N+1 : N;
      N       = alpha[violator] != 0.0 && aleph == 0.0 ? N-1 : N;
      E_alpha = E_alpha + fabs(aleph) - fabs(alpha[violator]);
   }

   alpha[violator] = aleph;
   set[violator]   = aleph == 0.0 ? I_z : I_nz;
}

/******************************************************************************

Procedure   : etherial

Parameters  :

Returns     :

Description :

******************************************************************************/

void etherial(const char *format, ...)
{
   char str[81];

   int i;

   va_list ap;

   va_start(ap, format);

   vsnprintf(str, 81, format, ap);

   mexPrintf(str);

   for (i = 0; i < 81 && str[i] != '\0'; i++)
   {
      str[i] = '\b';
   }

   mexPrintf(str);

   va_end(ap);
}

/******************************************************************************

Procedure   :

Parameters  :

Returns     :

Description :

******************************************************************************/


int i;

void blogreg()
{
   int i, violator;

   const char *fmt = "epoch = %-3d : Q = %-10.6f : E_w = %-10.6f : lambda = %f : N = %d  ";

   /* initialisation */

   epoch   = 0;
   E_alpha = 0.0;
   N       = 0;

   compute_objective();
  
   etherial(fmt, epoch++, Q, E_alpha, lambda, N);

   /* optimise model one parameter at a time */

   for (i = 0; i < 10000; i++)
   {
      /* find maximally violating inactive parameter */

      violator = findMaxViolator(I_nz); 

      if (violator != -1)
      {
         optimiseAlpha(violator);
      }
      else
      {
	 /* find maximally violating inactive parameter */

         violator = findMaxViolator(I_z); 

         if (violator != -1)
         {
            optimiseAlpha(violator);
         }

	 /* report progress */

	 compute_objective();

	 etherial(fmt, epoch++, Q, E_alpha, lambda, N);
      }

      /* update lambda */

      lambda = E_alpha == 0.0 ? lambda : N/E_alpha;

      if (violator == -1)
      {
         break;
      }
   }

   etherial("                                                                              ");
}

/******************************************************************************

Procedure   : mexFunction

Parameters  : int      nlhs   - number of output arguments
              mxArray *plhs[] - array of output arguments
	      int      nrhs   - number of input arguments
	      mxArray *prhs[] - array of input arguments

Returns     : void

Description : MEX gateway handling transfer of data to and from MATLAB, the
              calling sequence is something like

                 ALPHA = SLOGREG(X, T, LAMBDA, TOL, ALPHA)

              where X is the design matrix, T is a vector of target values,
	      GAMMA is the regularisation parameter, and ALPHA is a vector
	      of coefficients (this is optional as an input argument).

******************************************************************************/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   int i;

   double *X, *T;

   const mxArray *aleph = NULL;

   /* check number of input arguments */

   if (nrhs < 2 || nrhs > 4)
   {
      mexErrMsgTxt("wrong number of input arguments");
   }

   /* check number of output arguments */

   if (nlhs > 1)
   {
      mexErrMsgTxt("wrong number of output arguments");
   }

   /* get input patterns */

   ntp = mxGetM(prhs[0]); 
   d   = mxGetN(prhs[0]); 
   X   = (double*)mxGetPr(prhs[0]);

   /* get target patterns */

   if (mxGetM(prhs[1]) != ntp)
   {
      mexErrMsgTxt("T must have the same number of rows as X");
   }

   if (mxGetN(prhs[1]) != 1)
   {
      mexErrMsgTxt("T must be a column vector");
   }

   T = (double*)mxGetPr(prhs[1]);

   /* parse optional parameters */

   if (nrhs >= 3)
   {
      if (mxGetM(prhs[2]) == 1 && mxGetN(prhs[2]))
      {
         /* interpret as a value for tol */

         tol = mxGetScalar(prhs[2]);

	 if (tol < 0)
	 {
            mexErrMsgTxt("TOL must be a positive scalar");
	 }
      }
      else
      {
	 tol   = 1e-6;
         aleph = prhs[2];         
      }
   }

   if (nrhs == 4)
   {
      aleph = prhs[3];         
   }

   /* allocate matrix for result */

   plhs[0] = mxCreateDoubleMatrix(d, 1, mxREAL);
   alpha   = mxGetPr(plhs[0]);

   /* initialise alpha vector */

   if (aleph != NULL)
   {
      /* initial values supplied by user */

      if (mxGetN(aleph) != 1)
      {
         mexErrMsgTxt("ALPHA must be a column vector");
      }

      if (mxGetM(aleph) != d)
      {
         mexErrMsgTxt("ALPHA must have the same number of rows as X has columns");
      }

      memcpy(alpha, mxGetPr(aleph), d*sizeof(double));
   }
   else
   {
      /* initialise all alphas to zero */

      for (i = 0; i < d; i++)
      {
         alpha[i] = 0.0;
      }
   }

   /* initialisation */

   setvbuf(stdout, NULL, _IONBF, 0L);

   set    = (int*)     mxCalloc(d,   sizeof(int));
   xi     = (double*)  mxCalloc(ntp, sizeof(double));
   exp_xi = (double*)  mxCalloc(ntp, sizeof(double));
   F      = (double*)  mxCalloc(d,   sizeof(double));
   delta  = (double*)  mxCalloc(d,   sizeof(double));
   tmp1   = (double*)  mxCalloc(ntp, sizeof(double));
   tmp2   = (double*)  mxCalloc(ntp, sizeof(double));
   x      = (double**) mxCalloc(d,   sizeof(double*));
   t      = (double*)  mxCalloc(ntp, sizeof(double));

   for (i=0; i < d; i++)
   {
      x[i]     = &X[i*ntp];
      set[i]   = I_z;
      alpha[i] = 0.0;
   }

   set[0] = I_nz;

   for (i = 0; i < ntp; i++)
   {
      xi[i]     = 0.0;
      exp_xi[i] = 1.0;
      t[i]      = 2*T[i] - 1;
   }

   /* get on with it */

   blogreg();

   /* bye bye... */
}

/***************************** That's all Folks! *****************************/

