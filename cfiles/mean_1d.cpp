#include "mex.h"
#include "math.h"
/* Compute the arithmetic mean rapidly */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

int ii;
double media;
double *x;
int dimx;
double *out;

 dimx=mxGetNumberOfElements(prhs[0]);
 media=0;
 x=(double *)mxGetPr(prhs[0]);
 for (ii=0;ii<dimx;ii++)
  { media=media+*(x+ii);
  }
 media=media/dimx;

plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
out=mxGetPr(plhs[0]);
*(out)=media;
}