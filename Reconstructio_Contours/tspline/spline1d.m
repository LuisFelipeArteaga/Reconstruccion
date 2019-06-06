function [y, l, A] = spline1d (x_out, x_data, y_data, x_slope, y_slope, t, cutoff)
%	$Id: spline1d.m,v 1.2 2014/06/29 03:31:27 myself Exp $
% SPLINE1D	1-D interpolation using Green's function for a spline in tension
%
%	SPLINE1D will find a spline-based curve using continuous curvature splines
%	in tension (if set).  The algorithm uses the Green's function for the spline.
%	You can supply data constrains, slope constrains, or a mix of both.
%	Solution can be evaluated at arbitrary locations
%
%	Typically you use one of the following 3 call formats:
%
%	y = spline1d (x_out, x_data, y_data, x_slope, y_slope)
%	y = spline1d (x_out, x_data, y_data, x_slope, y_slope, t)
%	y = spline1d (x_out, x_data, y_data, x_slope, y_slope, t, cutoff)
%
% The input parameters are:
%
%	x_out	- Desired output x positions
%
%	x_data	- coordinates of points with data constraints
%	y_data	- data constraints at the above points
%	x_slope	- coordinates of points with slope constraints
%	y_slope	- slope constraints at the above points
%
%	t	- tension to use, 0 <= t <= 1
%		  if t is a vector of length 2 the second value is taken as the lengthscale
%	cutoff	- if set, eigenvalues whose ratio to the maximum eigenvalue are smaller than
%		  cutoff are zeroed out before the curve is evaluated.
%
%	One of (x_data, y_data) and (x_slope, y_slope) can be ([], [])
%	t, if not set, defaults to 0 (cubic spline).  t = 1 gives linear interpolation
%
% The output values are:
%
%	y	- the interpolation
%	l	- optionally, the eigenvalues of the linear system
%
%  See Wessel, P, D. Bercovici, 1998, Gridding with Splines in Tension : A Green function Approach,
%	Math. Geol., 30, 77-93.

% Pick a reasonable(?) lengthscale

length_scale = (max(x_out) - min(x_out)) / 50;

if (nargin == 5)	% No tension selected, set default
	t = 0;
else
	nt = length(t);
	if (nt == 2)	% User gave both tension and lengthscale
		length_scale = t(2);
		t = t(1);
	end
	if (t < 0.0 | t > 1.0)
		error ('spline1d: tension must be 0 <= t <= 1 !')
	end
end
if (nargin < 7)	% cutoff not set
	cutoff = 0;
end

% Misc initializations

if (t < 1)
	p = sqrt (t / (1 - t));
	p = p / length_scale;
end
n0 = 0;
n1 = 0;

% First we must enforce the use of column vectors for the data constrainst

[m,n] = size (x_out); if (m < n), x_out = x_out'; end
[m,n] = size (x_data); if (m < n), x_data = x_data'; end
[m,n] = size (y_data); if (m < n), y_data = y_data'; end
[n0,m0] = size (y_data);
[m,n] = size (x_slope); if (m < n), x_slope = x_slope'; end
[m,n] = size (y_slope); if (m < n), y_slope = y_slope'; end
n1 = length (x_slope);

% Assembly final xp, yp vectors (possibly combination of data and slopes)

xp = [x_data; x_slope];
yp = [y_data; y_slope];

% Now build the square n x n linear system that must be solved for the alpha's

n = n0 + n1;
A = zeros (n, n);
for i = 1:n0		% First add equations for data constraints
	r = xp(i) - xp;
	ar = abs(r);
	if (t == 0)
		A(i,:) = (ar .^ 3)';
	elseif (t == 1)
		A(i,:) = ar';
	else
		A(i,:) = (exp(-p * ar) + p * ar)';
	end
end
for i = 1:n1		% Then add equations for slope constraints
	j = i + n0;
	r = xp(j) - xp;
	ar = abs(r);
	if (t == 0)
		A(j,:) = 3.0 * (r .* ar)';
	elseif (t == 1)
		A(j,:) = sign (r)';
	else
		A(j,:) = p * (1.0 - exp(-p * ar))';
	end
end

% Done building square linear system, now solve it

if (cutoff > 0.0)	% Solve using SVD
	[U, S, V] = svd (A);
	s = diag (S);
	if (nargout == 2)	% Return eigenvalues
		l = s;
	end
	k = find ((s / s(1)) < cutoff);
	s = s .^ (-1);
	s(k) = zeros (size(k));
	S = diag (s);
	alpha = (V * S * U') * yp;
else			% Solve directly
	alpha = A \ yp;
	if (nargout == 2)	% Return eigenvalues
		l = svd(A);
	end
end

% Now evaluate final solution at output locations

y = zeros (length(x_out),m0);
for i = 1:n
	r = xp(i) - x_out;
	ar = abs(r);
	if (t == 0)
		y = y + (ar .^ 3) * alpha(i,:);
	elseif (t == 1)
		y = y + ar * alpha(i,:);
	else
		y = y +  (exp(-p * ar) + p * ar) * alpha(i,:);
	end
end

