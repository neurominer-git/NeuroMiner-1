#include "mex.h"
#include "math.h"
#include "matrix.h"
//#include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
int errfl = 0;
double *Y = mxGetPr(prhs[0]);
double *iX = mxGetPr(prhs[1]);
double *out;
int mY = mxGetM(prhs[0]);
int nY = mxGetN(prhs[0]);
int miX = mxGetM(prhs[1]);
int niX = mxGetN(prhs[1]);
//int numthread = omp_get_num_threads();
int numthread = mxGetScalar(prhs[2]);
int chunk = nY/numthread;
if (nrhs != 3){
  errfl=1;
}
if (niX > 1){
  errfl = 1;
}
// if (mY != miY)
//   errfl = 1;
if (miX > mY){
  errfl = 1;
}

if (errfl == 1) {
  mexPrintf("\nRows of Y = %d; Cols of Y = %d",mY,nY);
  mexPrintf("\nRows of iX = %d; Cols of iX = %d\n",miX,niX);
  mexErrMsgTxt("Usage: extrY(<m x n data matrix>, <p x 1> index of rows to extract");
}

//Output matrix
plhs[0]=mxCreateDoubleMatrix(miX,nY,mxREAL);
out=mxGetPr(plhs[0]);

#pragma omp parallel private shared(iX, Y, out, nY, mY, miX)
{
  int i, j, ind;
  #pragma omp for schedule(static, chunk)
  for (i = 0; i < nY; ++i)
  {
    for (j = 0; j < miX; ++j)
    {
      ind = (int)iX[j]-1;
      out[i*miX+j] = Y[i*mY+ind];
    }
  }
}
}