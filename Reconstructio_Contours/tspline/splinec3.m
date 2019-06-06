function [x, y, z, l] = splinec3 (xi, yi, zi, data, t, cutoff, n_out)
%	$Id: splinec3.m,v 1.1.1.1 2005/01/26 00:45:57 pwessel Exp $
% SPLINEC3	Curve-in-space interpolation using Green's function for a spline in tension
%
%	SPLINEC3 will find a spline-based curve using continuous curvature splines
%	in tension (if set).  The algorithm uses the Green's function for the spline.
%	You can supply data constrains, slope constrains, or a mix of both.
%	Solution is evaluated at equidistant locations along the resulting curve in
%	the x-y-z domain
%
%	[x, y, z]    = splinec3 (xi, yi, zi, data, t, cutoff, n_out)
%	[x, y, z, l] = splinec (xi, yi, zi, data, t, cutoff, n_out)
%
% The input parameters are:
%
%	xi, yi, zi	- coordinates of data constraints
%	data	- 9999 for regular points and angle in 0-360 degrees for slopes
%		 if all points are x,y only, data can be a constant 9999
%
%	t	- tension to use, 0 <= t <= 1
%		  if t is a vector of length 2 the second value is taken as the lengthscale
%	cutoff	- if set, eigenvalues whose ratio to the maximum eigenvalue are smaller than
%		  cutoff are zero'ed out before the curve is evaluated.
%
%
%	n_out	- Number of equidistant output points along the curve.  These points
%		  are augmented with the original data locations to arrive at the final
%		  output points; done to ensure the points will go through constraints
%
%	t, if not set, defaults to 0 (cubic spline).  t = 1 gives linear interpolation
%
% The output values are:
%
%	x, y, z	- the interpolating spatial curve
%	l	- optionally, the eigenvalues of the linear system
%
%  See Wessel, P, D. Bercovici, 1998, Gridding with Splines in Tension : A Green function Approach,
%       Math. Geol., 30, 77-93.


% First we must enforce the use of column vectors for the data constrainst

[m,n] = size (xi); if (m < n), xi = xi'; end
[m,n] = size (yi); if (m < n), yi = yi'; end
[m,n] = size (zi); if (m < n), zi = zi'; end
[m,n] = size (data); if (m < n), data = data'; end
n = length (yi);

if (length(data) == 1)	% All points
	k_pts = 1:n;
	k_slp = [];
else
	k_pts = find (data == 9999);
	k_slp = find (data ~= 9999);
end
data = data * pi / 180;

% Calculate distances

s = zeros (size (xi));
for i=2:n
	s(i) = s(i-1) + sqrt ((xi(i) - xi(i-1))^2 + (yi(i) - yi(i-1))^2 + (zi(i) - zi(i-1))^2);
end

% Now evaluate final solution at output locations that include the
% input locations as well as n_out new locations

s_out = sort ([linspace(0,s(n),n_out)'; s]);
ds = [1; diff(s_out)];
use = find (ds > 0);
s_out = s_out(use);

% make the call to spline1d with 3 rhs

if (nargout == 3)
	xyz = spline1d (s_out, s(k_pts), [xi(k_pts) yi(k_pts) zi(k_pts)], s(k_slp), [cos(data(k_slp)) sin(data(k_slp))], t, cutoff);
elseif (nargout == 4)
	[xyz, l] = spline1d (s_out, s(k_pts), [xi(k_pts) yi(k_pts) zi(k_pts)], s(k_slp), [cos(data(k_slp)) sin(data(k_slp))], t, cutoff);
else
	[xyz, l, A] = spline1d (s_out, s(k_pts), [xi(k_pts) yi(k_pts) zi(k_pts)], s(k_slp), [cos(data(k_slp)) sin(data(k_slp))], t, cutoff);
end
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
