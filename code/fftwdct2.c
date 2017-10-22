/* compile using
mex libfftw3.a fftwdct2.c
*/

#include "mex.h"
#include <fftw3.h>
#include <math.h>
#include "matrix.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs,  const mxArray *prhs[])
{
    if( nrhs != 3)
    {
        exit(1);
    }
    unsigned int m,n;
    m = (unsigned int)mxGetScalar(prhs[0]);
    n = (unsigned int)mxGetScalar(prhs[1]);
    
    
    unsigned int N = m*n;
    
    double *in = (double *)mxGetPr(prhs[2]);
    
    plhs[0] =  mxCreateDoubleMatrix(N,1,mxREAL);
    double *out = (double *)mxGetPr(plhs[0]);
    
    
    
    fftw_plan p = fftw_plan_r2r_2d(n, m,in,out, FFTW_REDFT10, FFTW_REDFT10,
                                   FFTW_ESTIMATE);
    
    fftw_execute(p);
    fftw_destroy_plan(p);
    
    /* scale the matrix out */
    double scale1 = 0.25/sqrt(N);
    out[0] *= scale1;
    
    scale1 *= 2.0; /* 1/(2*sqrt(m*n)) */
    double scale2 = scale1 * M_SQRT1_2;
    int i,j;
    
    
    for(i=1;i<n;i++)
    {
        out[i*m] *= scale2;
        for(j=1;j<m;j++)
            out[i*m + j] *= scale1;
    }
    
    
    
    for(i=1;i<m;i++) 
        out[i] *= scale2; /* first row */

    return ;
}

