#include "mex.h"
#include "math.h"
#include "matrix.h"

inline double intpow( double base, int exponent )
{
double curpwr = base;
int i;
for(i=1; i<=exponent; ++i)
	curpwr = base * curpwr;
return curpwr;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
double *out, *in;
int i,j, M, N;
in = mxGetPr(prhs[0]);
M = mxGetM(prhs[0]);
N = mxGetN(prhs[0]);

for (i=0;i<N;++i)
{
	for (j=0;j<M;++j)
		mexPrintf("\nRow %d, Col %d: Val=%g",j,i,in[i*M+j]);
}

plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
out = mxGetPr(plhs[0]);

}