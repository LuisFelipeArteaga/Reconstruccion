#	$Id: README.TXT,v 1.2 2014/06/29 03:34:21 myself Exp $
# README file for tspline distribution

This distribution contains the Matlab (and optional C mex subfunctions) developed
for the paper

Wessel, P., and D. Bercovici, 1998, Gridding with Splines in Tension : A Green
   function Approach, Math. Geol., 30, 77-93.

Please cite the paper if you use these functions in your published research.

The files supplied are:

spline1d.m		1-D spline in tension
spline2d.m		2-D spline in tension
spline3d.m		3-D spline in tension
splinec2.m		2-D curve-in-plane spline in tension
splinec3.m		3-D curve-in-space spline in tension
------------------------------------------------------------
spline2dgreen.m		Sub-function used by spline2d.m
spline2dgrad.m		Sub-function used by spline2d.m
detrend2.m		Sub-function used by spline2d.m
spline2d_grad.c		C mex function that can replace spline2dgrad.m
spline2d_grad.m		Dummy m file for documentation when C mex is used
spline2d_green.c	C mex function that can replace spline2dgreen.m
spline2d_green.m	Dummy m file for documentation when C mex is used
Makefile		Old Makefile showing how to build mex

Notes: 
1. All files are set to use the distributed Matlab subfunctions.  For additional
   speedup you can compile the two *.c files via Matlab's mex building machinery.
   The Makefile I supply is very old and you are likely to have to figure out how
   to build the mex files for your system.  If you do this successfully you will
   also need to edit spline2d.m and have it call spline2d_green (the C mex version)
   instead of spline2dgreen (the Matlab version); same for spline2d_grad.c
2. All the above splines as well as additional splines are available via the greenspline
   module in GMT, the Generic Mapping Tools (gmt.soest.hawaii.edu).

Paul Wessel
June 28, 2014
