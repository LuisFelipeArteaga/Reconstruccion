function [w, l] = spline3d (x_out, y_out, z_out, xi, yi, zi, wi, t)
%	$Id: spline3d.m,v 1.2 2014/06/29 03:31:27 myself Exp $
% SPLINE3D	3-D interpolation using Green's function for splines in tension
%
%	SPLINE3D will find a spline-based hyper-surface using continuous curvature splines
%	in tension (if set).  The algorithm uses the Green's function for the spline.
%	Solution can be evaluated on a 2-d horizontal plane using equidistant or arbitrary locations
%
%	Use the following call formats:
%
%	w = spline3d (x_out, y_out, z_out, xi, yi, zi, wi, t)
%
% The input parameters are:
%
%	x_out	- Desired output x positions (can be vector or matrix)
%	y_out	- Desired output y positions (can be vector or matrix)
%	z_out	- Desired output z level (a constant)
%
%	xi	- x-coordinates of points with data constraints
%	yi	- y-coordinates of points with data constraints
%	zi	- z-coordinates of points with data constraints
%	wi	- data constraints at the above points
%	t	- tension to use, 0 <= t < 1
%		  if t is a vector of length 2 the second value is taken as the lengthscale
%
% The output values are:
%
%	w	- the interpolation
%	l	- optionally, the eigenvalues of the linear system
%
%  See Wessel, P, D. Bercovici, 1998, Gridding with Splines in Tension : A Green function Approach,
%       Math. Geol., 30, 77-93.

% Pick a reasonable(?) lengthscale

IM = sqrt (-1);
length_scale = sqrt ((max(x_out(:)) - min(x_out(:)))^2 + (max(y_out(:)) - min(y_out(:)))^2) / 50;


% Misc initializations

p = sqrt (t / (1 - t));
p = p / length_scale;
if (t == 0)
	p = 1;
end

% First we must enforce the use of column vectors for the data constrainst

[m,n] = size (xi); if (m < n), xi = xi'; end
[m,n] = size (yi); if (m < n), yi = yi'; end
[m,n] = size (zi); if (m < n), zi = zi'; end
[m,n] = size (wi); if (m < n), wi = wi'; end
n = length (wi);

% Now build the square n x n linear system that must be solved for the alpha's

%disp ('build matrix')
%tic
A = zeros (n, n);
for i = 1:n
	u = p * sqrt ((xi(i) - xi) .^2 + (yi(i) - yi) .^2 + (zi(i) - zi) .^2);
	k = find (u > 0.0);
	A(i,:) = -ones(1,n);
	A(i,k) = ((exp (-u(k)) - 1.0) ./ u(k))';
%	A(i,:) = u';
end
%toc
%disp ('solve matrix')
%tic
% Done building square linear system, now solve it

alpha = A \ wi;
%toc
if (nargout == 2)	% Return eigenvalues
%	disp ('find eigenvalues')
%	tic
	l = svd (A);
%	toc
end

%disp ('evaluate')
%tic
% Now evaluate final solution at output locations

w = zeros (length(y_out),length(x_out));
[X, Y] = meshgrid(x_out, y_out);
for i = 1:length(alpha)
	u = p * sqrt ((X(:) - xi(i)) .^2 + (Y(:) - yi(i)) .^2 + (z_out - zi(i)) .^2)+eps;
	k = find (u > 0.0);
	v = -alpha(i) * ones(size(u));
	v(k) = alpha(i) * ((exp (-u(k)) - 1.0) ./ u(k));
%	v = alpha(i) * u;
	w(:) = w(:) + v;
end

%toc
%disp ('done')
