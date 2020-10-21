#include "mex.h"
#include "math.h"
#include "matrix.h"
//#include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
int errfl = 0;
double *Y = mxGetPr(prhs[0]);
int dim = (int)mxGetScalar(prhs[1]);
double *min;
double *minind;
int mY = mxGetM(prhs[0]);
int nY = mxGetN(prhs[0]);
//int numthread = omp_get_num_threads();
int numthread = mxGetScalar(prhs[2]);
int chunk;
int i, j;
double minY, minI;

if (nrhs != 3){
  errfl=1;
}

if ( dim == 1){ // across rows
    //Output matrix
    plhs[0]=mxCreateDoubleMatrix(1,nY,mxREAL);
    min = mxGetPr(plhs[0]);
    plhs[1]=mxCreateDoubleMatrix(1,nY,mxREAL);
    minind = mxGetPr(plhs[1]);
    
    #pragma omp parallel shared(Y, min, minind, nY, mY)
    {
    chunk = nY/numthread;
	#pragma omp for schedule(static, chunk)
	for (i = 0; i < nY; ++i)
	{
	      minY = Y[i*mY]; minI = 0;
	      //mexPrintf("\nY=%g",Y[i*Jx]);
	      for (j = 1; j < mY; ++j)
	      {
		  //mexPrintf("\nY=%g",Y[i*Jx+j]);
		  if ( Y[i*mY+j] < minY ){
		      minY = Y[i*mY+j];
		      minI = j;
		  }
	      }
	      min[i] = minY;
	      minind[i] = (double)minI + 1 ;
	}
    }

}
else if ( dim == 2 ){ // across columns

    //Output matrix
    plhs[0]=mxCreateDoubleMatrix(mY,1,mxREAL);
    min = mxGetPr(plhs[0]);
    plhs[1]=mxCreateDoubleMatrix(mY,1,mxREAL);
    minind = mxGetPr(plhs[1]);
    #pragma omp parallel shared(Y, min, minind, nY, mY)
    {
        chunk = mY/numthread;
	#pragma omp for schedule(static, chunk)
	for (i = 0; i < mY; ++i)
	{
	      minY = Y[i]; minI = 0;
	      //mexPrintf("\nY=%g",Y[i*Jx]);
	      for (j = 1; j < nY; ++j)
	      {
		  //mexPrintf("\nY=%g",Y[i*Jx+j]);
		  if ( Y[j*mY+i] < minY ){
		      minY = Y[j*mY+i];
		      minI = j;
		  }
	      }
	      min[i] = minY;
	      minind[i] = (double)minI + 1;
	}
    }
}
else {
    mexErrMsgTxt("You can either find the min across rows (flag = 1) or across columns (flag = 2)");
}
}