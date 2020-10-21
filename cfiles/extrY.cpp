#include "mex.h"
#include "math.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
int errfl = 0;
double *Y = mxGetPr(prhs[0]);
double *iY = mxGetPr(prhs[1]);
double *out;
int mY = mxGetM(prhs[0]);
int nY = mxGetN(prhs[0]);
int miY = mxGetM(prhs[1]);
int niY = mxGetN(prhs[1]);
int chunk = niY/4;
if (nrhs > 2){
  errfl=1;
}
if (miY > 1){
  errfl = 1;
}
// if (mY != miY)
//   errfl = 1;
if (niY > nY){
  errfl = 1;
}

if (errfl == 1) {
  mexPrintf("\nRows of Y = %d; Cols of Y = %d",mY,nY);
  mexPrintf("\nRows of iY = %d; Cols of iY = %d\n",miY,niY);
  mexErrMsgTxt("Usage: extrY(<m x n data matrix>, <1 x n> index of columns to extract");
}

//Output matrix
plhs[0]=mxCreateDoubleMatrix(mY,niY,mxREAL);
out=mxGetPr(plhs[0]);

#pragma omp parallel shared(Y, out, mY, niY)
{
  int i, j, ind;
  #pragma omp for schedule(static, chunk)
  for (i = 0; i < niY; ++i)
  {
    for (j = 0; j < mY; ++j)
    {
      ind = iY[i]-1;
      out[i*mY+j] = Y[ind*mY+j];
    }
  }
}
}