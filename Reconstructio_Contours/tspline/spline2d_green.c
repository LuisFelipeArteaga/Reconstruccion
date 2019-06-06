/*
 * Program:	spline2d_green.c
 * Purpose:	matlab callable routine to calculate Green function for 2-d
 *		spline in tension
 * Author:	P Wessel
 */
 
#include <math.h>
#include "/usr/local/matlab4/extern/include/mex.h"

 /* Gateway routine */
   
void mexFunction (int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[])
{
	double *x, *cc, *k0, c;
	int n, m, np;
	void spline2d_green();
 
	if (nrhs != 2 || nlhs > 1) {
		mexPrintf (" usage: k0 = spline2d_green (x, c);\n");
		return;
	}

	m = mxGetM (prhs[0]);
	n = mxGetN (prhs[0]);
	np = n * m;

	/* Get pointers to input arrays */

	x = mxGetPr (prhs[0]);
	cc = mxGetPr (prhs[1]);
	c = cc[0];

	/* Create a matrix for the return array */

	plhs[0] = mxCreateFull (m, n, REAL);
    
	k0 = mxGetPr (plhs[0]);
 
	/* Do the actual computations in a subroutine */
 
	spline2d_green (x, c, k0, np);
}


/* spline2d_green computes the Green function for a 2-d spline possibly
 * in tension, G(u) = G(u) - log(u), where u = c * x and c = sqrt (t/(1-t)).
 * The modified Bessel function K of order zero is based on Num. Rec.
 * All x must be >= 0.  When c = 0 it degenerates to x^2 * log(x)
 */

void spline2d_green (x, c, G, n)
double x[], c, G[];
int n;
{
	int i;
	double y, z, ic, cx, g0, t;

	if (c == 0.0) {	/* Just regular spline */
		for (i = 0; i < n; i++) G[i] = (x[i] == 0.0) ? 0.0 : x[i] * x[i] * log (x[i]);
		return;
	}

	ic = 1.0 / c;
	g0 = 0.11593158056;	/* log(2) - 0.5772156 */

	for (i = 0; i < n; i++) {
		if (x[i] == 0.0) {
			G[i] = 0.0;
		}
		else if (x[i] <= 2.0 * ic) {
			cx = c * x[i];
			y = 0.25 * (t = cx * cx);
			z = t / 14.0625;
			G[i] = (-log(0.5*cx) * (z * (3.5156229 + z * (3.0899424 + z * (1.2067492 + z * (0.2659732 + z * (0.360768e-1 + z * 0.45813e-2))))))) + (y * (0.42278420 + y * (0.23069756 + y * (0.3488590e-1 + y * (0.262698e-2 + y * (0.10750e-3 + y * 0.74e-5))))));
		}
		else {
			y = 2.0 * ic / x[i];
			cx = c * x[i];
			G[i] = (exp (-cx) / sqrt (cx)) * (1.25331414 + y * (-0.7832358e-1 + y * (0.2189568e-1 + y * (-0.1062446e-1 + y * (0.587872e-2 + y * (-0.251540e-2 + y * 0.53208e-3)))))) + log (cx) - g0;
		}
	}
}
