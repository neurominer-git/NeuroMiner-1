#include "mex.h"
#include "math.h"
/* Compute the arithmetic mean rapidly */
double mean1d(double *x, int dimx);
double cov1d(double *x, int dimx, double avg);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
int ii;
int dimx1, dimx2;
double *x1; 
double *x2;
double *out;
double avg1, avg2, avgdiff, var1, var2, fdr;

dimx1=mxGetNumberOfElements(prhs[0]);
dimx2=mxGetNumberOfElements(prhs[1]);

x1=(double *)mxGetPr(prhs[0]);
x2=(double *)mxGetPr(prhs[1]);

/* Compute sums and means */
avg1 = mean1d(x1, dimx1);
var1 = cov1d(x1, dimx1, avg1);

avg2 = mean1d(x2, dimx2);
var2 = cov1d(x2, dimx2, avg2);
avgdiff = avg1-avg2;

fdr = (avgdiff*avgdiff)/(var1+var2);

plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
out=mxGetPr(plhs[0]);
*(out)=fdr;
}

double mean1d(double *x, int dimx){
int ii;
double media=0;
for (ii=0;ii<dimx;ii++)
  { media=media+*(x+ii);
  }
return media/dimx;
}

double cov1d(double *x, int dimx, double avg){
int ii;
double cov = 0;
double val;
for(ii=0;ii<dimx;ii++)
{
	val=*(x+ii)-avg;
	cov=cov+val*val;
}
return cov/(dimx-1);
}