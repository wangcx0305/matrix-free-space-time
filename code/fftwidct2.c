/* compile using
 * mex libfftw3.a fftwidct2.c
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
    int i,j;
    
    
    m = (unsigned int)mxGetScalar(prhs[0]);
    n = (unsigned int)mxGetScalar(prhs[1]);
    
    
    unsigned int N = m*n;
    
    double *in = (double *)mxGetPr(prhs[2]);
    
    plhs[0] =  mxCreateDoubleMatrix(N,1,mxREAL);

    double *out = (double *)mxGetPr(plhs[0]);
    
    
    double scale = 1/sqrt(N);
    

    out[0] = scale*in[0];
    
    scale *= 0.5;
    for(j=1;j<n;j++)
	for(i=1;i<m;i++)
            out[i + m*j] = scale*in[i + m*j];
    
    
    scale *= 1.41421356237309504880;

    for(i=1;i<m;i++)
        out[i] = scale*in[i]; /* first column */
    for(i=1;i<n;i++)   
        out[i*m] = scale*in[i*m]; /* first row */


    fftw_plan p = fftw_plan_r2r_2d(n, m,out,out, FFTW_REDFT01, FFTW_REDFT01,
            FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    
    return ;
}


