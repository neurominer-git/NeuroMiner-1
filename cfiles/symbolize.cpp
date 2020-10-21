#include "mex.h"
#include "math.h"
#include "matrix.h"

double intpow( double base, int exponent ){
if (exponent == 0)
	return 1;
else
{
	double curpwr = base;
	int i;
	for(i=1; i<exponent; ++i)
		curpwr = base * curpwr;
	return curpwr;
}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
// Input arguments: plhs[0] = Fea, plhs[1]= ns, plhs[2]= ne, plhs[3]= seql, plhs[4] = winx 
//

int M = mxGetM(prhs[0]);
int N = mxGetN(prhs[0]);
int ns = mxGetScalar(prhs[1]);
int ne = mxGetScalar(prhs[2]);
double seql = mxGetScalar(prhs[3]);
int winx = mxGetScalar(prhs[4]);
int Mseql = M-seql;
int ndiff = ne-ns+1;
double *binFea, *dFea, *Entropy;
double *Fea = mxGetPr(prhs[0]);

// symbolized data
plhs[0] = mxCreateDoubleMatrix(M,N,mxREAL);
dFea = mxGetPr(plhs[0]);

// bins used for symbolization: 
// row 1 = start
// row 2 = step
// row 3 = stop
plhs[1] = mxCreateDoubleMatrix(3,N,mxREAL);
binFea = mxGetPr(plhs[1]);

// Entropies 
plhs[2]=mxCreateDoubleMatrix(ndiff,N,mxREAL);
Entropy=mxGetPr(plhs[2]);

int chunk = N/4;
//mexPrintf("Chunk: %d",chunk);

#pragma omp parallel shared(Fea,dFea,binFea,Entropy,M,N,ns,ne,seql,winx,Mseql,ndiff)
{
int i, j, k, l, u, bestj = 0, flag, cnt;
double avg, var, val, ibin, datmin, datmax, ulim, olim, stdwin;
double Cx, Ex, maxEx, pEx, P;
double* C = new double[Mseql];
double* tFea = new double[M*ndiff];
#pragma omp for schedule(static,chunk)
for (i=0;i<N;++i)
{
	// Determine average and std
	avg = 0; var=0;
	for (j=0; j<M; ++j)
		avg += Fea[i*M+j];

	avg /= M;
	for (j=0; j<M; ++j)
	{
		val = Fea[i*M+j] - avg;
		var += val*val;
	}
	if (avg==0 && var ==0)
		continue;
	var = sqrt(var/(M-1));
	stdwin = winx*var;
	datmin = avg - stdwin;
	datmax = avg + stdwin;
	maxEx = 0;
	cnt=0;
	for (j=ns; j<=ne; ++j)
	{
		// Step 1: Determine partitioning
		ibin = (datmax-datmin)/j;
		//mexPrintf("\ni: %d, j: %d, ibin: %g", i, j, ibin);
		// Step 2: Digitize data
		for (k=0; k<M; ++k)
		{
			if (Fea[i*M+k] <= datmin) // Value below digitization window
			{
				tFea[cnt*M+k] = 0;
				continue;
			}
			else if (Fea[i*M+k] >= datmax) // Value above digitization window
			{
				tFea[cnt*M+k] = j;
				continue;
			}
			else
			{
				for (l=0; l<j; ++l) // Value within digitization window => determine optimal bin size
				{
					ulim = datmin + ibin*l; olim = datmin + ibin*(l+1);
					//mexPrintf("\n%g ? >=%g && <%g: %d", Fea[k], ulim, olim, l);
					if (Fea[i*M+k] >= ulim && Fea[i*M+k] < olim)
					{
						tFea[cnt*M+k] = l;
						continue;
					}
				}
			}
		}
	
		// Step 3: generate symbolic sequences
		l = 0;
		//mexPrintf("\n\nGenerate symbolic sequences.");
		while (l+seql < M)
		{
			k=l; Cx = 0;
			//mexPrintf("\nC:");
			for (u = seql-1; u>=0; --u)
			{
				Cx += intpow(j,u)*tFea[cnt*M+k];
				//mexPrintf("\nj=%d, u=%d, j^u=%g, dFea=%g, Cx=%g",j,u,intpow(j,u),dFea[i*M+k],Cx);
				++k;
			}
			C[l] = Cx;
			//mexPrintf("\nC: %g", C[l]);
			++l;
		}
		
		// Step 4: calculate Shannon's entropy of current sequence
		Ex = 0;
		//mexPrintf("\n\nCompute entropy.");
		for (k=0; k<Mseql; ++k)
		{
			flag=0;
			// Compute occurences of unique value
			if (k != 0)
			{
				for (l=0; l<k; ++l)
				{
					//mexPrintf("\nk=%d,l=%d => %g=%g?",k,l,C[k],C[l]);
					if ( C[k] == C[l] )
					{
						flag=1;
						//mexPrintf("\t found.");
					}
				}
			}
			if (flag == 0)
			{
				P = 1;
				for (l=k+1; l<Mseql; ++l)
				{
					if ( C[l] == C[k] )
						++P;
				}
				//mexPrintf("\nPorig = %g",P);
				P = P / Mseql;
				Ex = Ex + P*log(P);
				//mexPrintf(",P/Mseql=%g,Ex=%g",P,Ex);
			}
		}
		Ex = (-1/seql)*Ex;
		Entropy[i*ndiff+cnt] = Ex;
		//mexPrintf("\nEntropy=%g",Ex);
		if ( Ex > maxEx )
			maxEx = Ex;
	++cnt;
	}
	// Determine maximum according to the 99%-from-max criterion
	// that may be to simplistic, but it is fast compared to 
	// more sophisticated knee point detection algorithms

	for (j=0;j<ndiff;++j)
	{
		pEx = Entropy[i*ndiff+j]/maxEx;
		if (pEx > 0.99)
		{
			bestj = j;
			break;
		}
	}	
	ibin = (datmax-datmin)/(bestj+ns);
	binFea[i*3] = datmin;
	binFea[i*3+1] = ibin;
	binFea[i*3+2] = datmax;
	for (k=0; k<M; ++k)
		dFea[i*M+k] = tFea[bestj*M+k];
}
delete [] C;
delete [] tFea;
}
}
