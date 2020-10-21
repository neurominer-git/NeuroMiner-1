/* compile with: mex sp_factor_ratio.c
 * usage: Z = sp_factor_ratio(X, W, H)
 * compute Z = X ./ (W*H) but in a sparse manner
 * X must be sparse, cannot be full
 */
/* Written by Zhirong Yang
 */
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  double *X, *W, *H, *Z;
  mwIndex *irs,*jcs,ind,i,j,k,row;
  mwIndex   starting_row_index, stopping_row_index, current_row_index;
  mwSize col;
  mwSize m,n,nzmax,nzrows,r;
  double xhat;
  
  if ((nlhs < 1) || (nrhs!=3))
    mexErrMsgTxt("Usage: Z = sp_factor_ratio(X, W, H)");
  
  X = mxGetPr(prhs[0]);
  W = mxGetPr(prhs[1]);
  H = mxGetPr(prhs[2]);
  
  irs = mxGetIr(prhs[0]);
  jcs = mxGetJc(prhs[0]);

  m  = mxGetM(prhs[0]);
  n  = mxGetN(prhs[0]);
  r  = mxGetN(prhs[1]);

  plhs[0] = mxDuplicateArray(prhs[0]);
  Z  = mxGetPr(plhs[0]);
  
  ind = 0;
  for (col=0; col<n; col++) {
      starting_row_index = jcs[col];
      stopping_row_index = jcs[col+1];
      if (starting_row_index == stopping_row_index)
          continue;
      else {
          for (current_row_index = starting_row_index;
          current_row_index < stopping_row_index;
          current_row_index++)  {
              row = irs[ind];
              xhat = 0;
              for (k=0;k<r;k++){
                  xhat += W[m*k+row] * H[r*col+k];
              }
              Z[ind] = X[ind] / xhat;
              ind++;
          }
      }
  }
  
}
