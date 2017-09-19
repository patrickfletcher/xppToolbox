
/* Based on MatLab's example "yprime"
*
* This is a template MEX function for computing a function of the form:
*     [yp, aux]=f(t, y, p)
* such as those required for numerically integrating ODEs.  t[scalar], y[vector], p[vector]. 
*     yp = y prime (the derivative of the vector field given by f, at point (t,y) parameterized by the parameter vector p.
*     aux = vector of auxiliary output functions defined in and returned by f (xppaut's "aux" declarations)
*
* This file will be appended to a section of generated code containing:
*
* #include <math.h>
* #include "mex.h"
* user defines
* static void yprime(double	*t, double y[], double p[], double dy[], double aux[]);
*
* (C) Patrick Fletcher,  October 2015
*/

//ONLY DOUBLE WORKS using this template...
		   
/* Input Arguments */
#define	T_IN	prhs[0]
#define	Y_IN	prhs[1]
#define	P_IN	prhs[2]

/* Output Arguments */
#define	YP_OUT	plhs[0]
#define	AUX_OUT	plhs[1]

/* Other useful macros */
#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#if !defined(heav)
#define heav(x) ((x) > 0.0 ? 1.0 : 0.0) 
#endif

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    double *t,*y,*p,*yp,*aux; 
	size_t m,n;
    
    /* Check for proper number of arguments */
    if (nrhs != 3) { 
	    mexErrMsgIdAndTxt( "MATLAB:yprime:invalidNumInputs",
                "Three input arguments required."); 
    }
	
	if (nlhs > 2) {
	    mexErrMsgIdAndTxt( "MATLAB:yprime:maxlhs",
                "Too many output arguments."); 
    }
	
    /* Check the dimensions of T.  T should be scalar */ 
	if (!mxIsScalar(T_IN)) { 
	    mexErrMsgIdAndTxt( "MATLAB:yprime:invalidT",
                "YPRIME requires that T be a scalar single or double"); 
    } 

    /* Check the dimensions of Y.  Y should be N_VAR X 1 */ 
    m = mxGetM(Y_IN); 
    n = mxGetN(Y_IN);
	if (MAX(m,n) != N_VAR || MIN(m,n) != 1) { 
	    mexErrMsgIdAndTxt( "MATLAB:yprime:invalidY",
                "YPRIME requires that Y be a %i x 1 vector of single or double.", N_VAR); 
    } 
	
    /* Check the dimensions of P.  P should be N_PAR X 1 */ 
    m = mxGetM(P_IN); 
    n = mxGetN(P_IN);
	if (N_PAR>0 && (MAX(m,n) != N_PAR || MIN(m,n) != 1)) { 
	    mexErrMsgIdAndTxt( "MATLAB:yprime:invalidP",
                "YPRIME requires that P be a %i x 1 vector of singles or doubles.", N_PAR); 
    } 


	/* Input Type checking */
	if (!mxIsSingle(T_IN) && !mxIsDouble(T_IN)) { 
	    mexErrMsgIdAndTxt( "MATLAB:yprime:invalidT",
                "YPRIME requires that T be a scalar single or double"); 
    } 
	if (!mxIsSingle(Y_IN) && !mxIsDouble(Y_IN)) { 
	    mexErrMsgIdAndTxt( "MATLAB:yprime:invalidY",
                "YPRIME requires that Y be a %i x 1 vector of singles or doubles.", N_VAR); 
    } 
	if (!mxIsSingle(P_IN) && !mxIsDouble(P_IN)) { 
	    mexErrMsgIdAndTxt( "MATLAB:yprime:invalidP",
                "YPRIME requires that P be a %i x 1 vector of singles or doubles.", N_PAR); 
    } 

    /* Create a matrix for the return argument */ 
    YP_OUT = mxCreateDoubleMatrix( (mwSize)N_VAR, (mwSize)1, mxREAL); 
    AUX_OUT = mxCreateDoubleMatrix( (mwSize)N_AUX, (mwSize)1, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    yp = mxGetPr(YP_OUT);
	aux = mxGetPr(AUX_OUT);
    
    t = mxGetPr(T_IN);
    y = mxGetPr(Y_IN);
    p = mxGetPr(P_IN);
        
    /* Do the actual computations in a subroutine */
    yprime(t,y,p,yp,aux);


    return;
    
}
