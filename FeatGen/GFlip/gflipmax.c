#include "mex.h"
#include "math.h"
#include "matrix.h"
//#include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
int errfl = 0;
double *Y = mxGetPr(prhs[0]);
int dim = (int)mxGetScalar(prhs[1]);
double *out;
int mY = mxGetM(prhs[0]);
int nY = mxGetN(prhs[0]);
//int numthread = omp_get_num_threads();
int numthread = mxGetScalar(prhs[2]);
int i, j;
int chunk;
double maxY;

if (nrhs != 3){
  errfl=1;
}
if ( dim == 1){ // across rows
    //Output matrix
    plhs[0]=mxCreateDoubleMatrix(1,nY,mxREAL);
    out=mxGetPr(plhs[0]);
    #pragma omp parallel shared(Y, out, nY, mY)
    {
	chunk = nY/numthread;
	#pragma omp for schedule(static, chunk)
	for (i = 0; i < nY; ++i)
	{
	      maxY = Y[i*mY];
	      //mexPrintf("\nY=%g",Y[i*Jx]);
	      for (j = 1; j < mY; ++j)
	      {
		  //mexPrintf("\nY=%g",Y[i*Jx+j]);
		  if ( Y[i*mY+j] > maxY )
		      maxY = Y[i*mY+j];
	      }
	      out[i] = maxY;
	}
    }

}
else if ( dim == 2 ){ // across columns

    //Output matrix
    plhs[0]=mxCreateDoubleMatrix(mY,1,mxREAL);
    out=mxGetPr(plhs[0]);
    #pragma omp parallel shared(Y, out, nY, mY)
    {
	chunk = mY/numthread;
	#pragma omp for schedule(static, chunk)
	for (i = 0; i < mY; ++i)
	{
	      maxY = Y[i];
	      //mexPrintf("\nY=%g",Y[i*Jx]);
	      for (j = 1; j < nY; ++j)
	      {
		  //mexPrintf("\nY=%g",Y[i*Jx+j]);
		  if ( Y[j*mY+i] < maxY )
		      maxY = Y[j*mY+i];
	      }
	      out[i] = maxY;
	}
    }
}
else {
    mexErrMsgTxt("You can either find the min across rows (flag = 1) or across columns (flag = 2)");
}
}