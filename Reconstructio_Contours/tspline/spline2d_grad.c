/*
 * Program:	spline2d_grad.c
 * Purpose:	matlab callable routine to calculate gradient of Green function for 2-d
 *		spline in tension
 * Author:	P Wessel
 */
 
#include <math.h>
#include "/usr/local/matlab4/extern/include/mex.h"

 /* Gateway routine */
   
void mexFunction (int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[])
{
	double *dx, *dy, *i, *j, *cc, *dzdn, c, ni, nj;
	int n, m, np;
	void spline2d_grad();
 
	if (nrhs != 5 || nlhs > 1) {
		mexPrintf (" usage: dzdn = spline2d_grad (dx, dy, i, j, c);\n");
		return;
	}

	m = mxGetM (prhs[0]);
	n = mxGetN (prhs[0]);
	np = n * m;

	/* Get pointers to input arrays */

	dx = mxGetPr (prhs[0]);
	dy = mxGetPr (prhs[1]);
	i = mxGetPr (prhs[2]);
	j = mxGetPr (prhs[3]);
	cc = mxGetPr (prhs[4]);
	c = cc[0];
	ni = i[0];
	nj = j[0];

	/* Create a matrix for the return array */

	plhs[0] = mxCreateFull (m, n, REAL);
    
	dzdn = mxGetPr (plhs[0]);
 
	/* Do the actual computations in a subroutine */
 
	spline2d_grad (dx, dy, ni, nj, c, dzdn, np);
}


/* spline2d_grad computes the radial derivative of Green function for a 2-d spline possibly
 * in tension, G'(u) = 1/u - K1(u), where u = c * x and c = sqrt (t/(1-t)).
 * The modified Bessel function K of order one is based on Num. Rec.
 * All x must be >= 0.  When c = 0 it degenerates to x (2*log(x)-1)
 */

void spline2d_grad (dx, dy, ni, nj, c, dgdn, n)
double dx[], dy[], ni, nj, c, dgdn[];
int n;
{
	int i;
	double y, z, r, dgdr, ic, cx, cx2, icx;

	if (c == 0.0) {	/* Just regular spline */

		for (i = 0; i < n; i++) {
			r = hypot (dx[i], dy[i]);
			if (r == 0) {
				dgdn[i] = 0.0;
			}
			else {
				dgdr = 2.0 * log (r) - 1.0;
				dgdn[i] = dgdr * (dx[i] * ni + dy[i] * nj);
			}
		}
			
		return;
	}

	ic = 1.0 / c;

	for (i = 0; i < n; i++) {
		r = hypot (dx[i], dy[i]);
		if (r == 0.0) {
			dgdn[i] = 0.0;
			continue;
		}

		if (r <= 2.0 * ic) {
			cx = c * r;
			cx2 = cx * cx;
			y = 0.25 * cx2;
			z = cx2 / 14.0625;
			dgdr = -((log(0.5*cx) * (cx * (0.5 + z * (0.87890594 + z * (0.51498869 + z * (0.15084934 + z * (0.2658733e-1 + z * (0.301532e-2 + z * 0.32411e-3)))))))) + (1.0/cx) * (y * (0.15443144 + y * (-0.67278579 + y * (-0.18156897 + y * (-0.1919402e-1 + y * (-0.110404e-2 + y * (-0.4686e-4))))))));
		}
		else {
			cx = c * r;
			icx = 1.0 / cx;
			y = 2.0 * icx;
			dgdr = icx - ((exp (-cx) / sqrt (cx)) * (1.25331414 + y * (0.23498619 + y * (-0.3655620e-1 + y * (0.1504268e-1 + y * (-0.780353e-2 + y * (0.325614e-2 + y * (-0.68245e-3))))))));
		}
		dgdn[i] = c * dgdr * (dx[i] * ni + dy[i] * nj) / r;
	}
}
