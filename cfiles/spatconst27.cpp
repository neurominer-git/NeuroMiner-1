#include "mex.h"
#include "math.h"
#include "matrix.h"
//#define V(i,j,k) v[i+(j+k*N)*M]
//#define C(i,j,k) c[i+(j+k*N)*M]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
double *v;
double *c;
const mwSize *dimarray;
int chunk, X, Y, Z;
int i, j, k, l, m ,n;

mwSize dimnum;

dimarray = mxGetDimensions(prhs[0]);
dimnum = mxGetNumberOfDimensions(prhs[0]);

//mexPrintf("\n# dimensions: %d",dimnum);
//mexPrintf("\nM=%d, N=%d",M,N);
//for (i=0;i<dimnum;++i)
  //mexPrintf("\nDim %d: Num of elements %d", i, dimarray[i]);

v = mxGetPr(prhs[0]);
plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),mxGetDimensions(prhs[0]),mxDOUBLE_CLASS, mxREAL);
c = (double *)mxGetData(plhs[0]);

X=dimarray[0];
Y=dimarray[1];
Z=dimarray[2];

chunk = dimarray[0]/4;
//mexPrintf("Chunk: %d",chunk);
// Loop through all 3 dimensions
#pragma omp parallel private(i,j,k,l,m,n) shared(v, c, dimarray,X,Y,Z)
{
	double avg, cov, val;
	#pragma omp for
	for (i=1;i<X-1;++i)
	{
		for (j=1;j<Y-1;++j)
		{
			for (k=1;k<Z-1;++k)
			{
				// Extract values from 27k cube
				avg = 0; cov=0;
				//mexPrintf("\nCurrent position: x=%d, y=%d, z=%d.\nV=",i,j,k);
				for (n=k-1;n<=k+1;++n)
				{
					for (l=i-1;l<=i+1;++l)
					{
						for (m=j-1;m<=j+1;++m)
						{
							//mexPrintf(" / %d %g ",l+m*M+n*M,v[l+m*M+n*M]);
							//mexPrintf(" %g",v[l+m*M+n*N*M]);
							avg += v[l+m*X+n*X*Y];
						}
					}
				}
				avg /= 27;
				for (n=k-1;n<=k+1;++n)
				{
					for (l=i-1;l<=i+1;++l)
					{
						for (m=j-1;m<=j+1;++m)
						{
							val = v[l+m*X+n*X*Y]-avg;
							cov += val*val;
						}
					}
				}
				c[i+j*X+k*X*Y]=cov/26;
				//mexPrintf("\nVar: %g",cov/26);
			}
		}
	}
	//mexPrintf("\nfinished.");
}
}