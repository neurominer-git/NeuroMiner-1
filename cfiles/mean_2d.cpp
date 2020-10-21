#include "mex.h"
#include "math.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
int errfl = 0;

if (nrhs > 2){
  errfl=1;
}

double *Y = mxGetPr(prhs[0]);
int dim = mxGetScalar(prhs[1]);

int mY = mxGetM(prhs[0]);
int nY = mxGetN(prhs[0]);

//int chunk = nY/4;

if (errfl == 1) {
/*  mexPrintf("\nRows of Y = %d; Cols of Y = %d",mY,nY);
  mexPrintf("\nRows of iY = %d; Cols of iY = %d\n",miY,niY);*/
  mexErrMsgTxt("Usage: mean_2d(<m x n data matrix>, dimension to use [=1/=2]");
}

if (dim == 1){ // as in MATLAB
  //Output matrix
  plhs[0]=mxCreateDoubleMatrix(1,nY,mxREAL);
  double *out=mxGetPr(plhs[0]);
  //#pragma omp parallel shared(Y, out, mY, nY)
  //{
    int i, j, avg;
  //  #pragma omp for schedule(static, chunk)
    for (i = 0; i < nY; ++i){
      avg = 0;
      for (j = 0; j < mY; ++j){
	avg = avg + Y[i*mY+j];
      }
      out[i] = avg / mY; 
    } 
  //}
}
else if (dim == 2){
  plhs[0]=mxCreateDoubleMatrix(mY,1,mxREAL);
  double *out=mxGetPr(plhs[0]);
 // #pragma omp parallel shared(Y, out, mY, nY)
 // {
    int i, j, avg;
 //   #pragma omp for schedule(static, chunk)
    for (i = 0; i < mY; ++i){
	avg = 0;
	for (j = 0; j < nY; ++j){
	  avg = avg + Y[i*mY+j];
	}
	out[i] = avg / mY; 
    }
  //}
}


}