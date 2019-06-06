function [z, l] = spline2d (x_out, y_out, arg_3, arg_4, arg_5, arg_6, arg_7, arg_8, arg_9, arg_10, arg_11)
%	$Id: spline2d.m,v 1.3 2014/06/29 03:31:27 myself Exp $
% SPLINE2D	Gridding using Green's function for splines in tension
%
%	SPLINE2D will find a spline-based surface using continuous curvature splines
%	in tension (if set).  The algorithm uses the Green's function for the spline.
%	You can supply data constrains, slope constrains, or a mix of both.
%	Solution can be evaluated on a grid or at arbitrary locations
%
%	Use one of the following 6 call formats:
%
%	z = spline2d (x_out, y_out, x_data, y_data, z_data)
%	z = spline2d (x_out, y_out, x_data, y_data, z_data, t)
%	z = spline2d (x_out, y_out, x_slope, y_slope, i_slope, j_slope, z_slope)
%	z = spline2d (x_out, y_out, x_slope, y_slope, i_slope, j_slope, z_slope, t)
%	z = spline2d (x_out, y_out, x_data, y_data, z_data, x_slope, y_slope, i_slope, j_slope, z_slope)
%	z = spline2d (x_out, y_out, x_data, y_data, z_data, x_slope, y_slope, i_slope, j_slope, z_slope, t)
%
% The input parameters are:
%
%	x_out	- Desired output x positions (can be vector or matrix)
%	y_out	- Desired output y positions (can be vector or matrix)
%
%	x_data	- x-coordinates of points with data constraints
%	y_data	- y-coordinates of points with data constraints
%	z_data	- data constraints at the above points
%	x_slope	- x-coordinates of points with slope constraints
%	y_slope	- y-coordinates of points with slope constraints
%	i_slope	- x-component of slope direction for the above points
%	j_slope	- y-component of slope direction for the above points
%	z_slope	- slope constraints at the above points
%
%	t	- tension to use, 0 <= t < 1
%		  if t is a vector of length 2 the second value is taken as the lengthscale
%
% The output values are:
%
%	z	- the interpolation
%	l	- optionally, the eigenvalues of the linear system
%
%  See Wessel, P, D. Bercovici, 1998, Gridding with Splines in Tension : A Green function Approach,
%       Math. Geol., 30, 77-93.

% By default we call the Matlab functions spline2dgrad and spline2dgreen; alternatively you can compile
% spline2d_grad.c and spline2d_green.c into mex files and add in those underscores below in the calls.

% Pick a reasonable(?) lengthscale
verbose = 1;

IM = sqrt (-1);
length_scale = abs (max(x_out(:)) - min(x_out(:)) + IM * (max(y_out(:)) - min(y_out(:)))) / 50;

if (nargin == 5 | nargin == 10 | nargin == 7)	% No tension selected, set default
	t = 0;
	n_args = nargin;
else
	n_args = nargin - 1;
	t = eval (['arg_' int2str(nargin)]);
	if (length(t) == 2)	% User gave both tension and lengthscale
		length_scale = t(2);
		t = t(1);
	end
	if (t < 0.0 | t >= 1.0)
		error ('spline2d: tension must be 0 <= t < 1 !')
	end
end

% Now figure out what was passed and assign the data accordingly

if (n_args == 5 | n_args == 10)	% z_data supplied
	x_data = arg_3;
	y_data = arg_4;
	z_data = arg_5;
end
if (n_args == 7)	% only z_slope supplied
	x_slope = arg_3;
	y_slope = arg_4;
	i_slope = arg_5;
	j_slope = arg_6;
	z_slope = arg_7;
elseif (n_args == 10)	% z_slope supplied
	x_slope = arg_6;
	y_slope = arg_7;
	i_slope = arg_8;
	j_slope = arg_9;
	z_slope = arg_10;
end

% Misc initializations

p = sqrt (t / (1 - t));
p = p / length_scale;
n0 = 0;
n1 = 0;

% First we must enforce the use of column vectors for the data constrainst

if (n_args == 5 | n_args == 10)	% z_data supplied; check if we must transpose
	[m,n] = size (x_data); if (m < n), x_data = x_data'; end
	[m,n] = size (y_data); if (m < n), y_data = y_data'; end
	[m,n] = size (z_data); if (m < n), z_data = z_data'; end
	n0 = length (z_data);
end
if (n_args == 7 | n_args == 10)	% z_slope supplied; check if we must transpose
	[m,n] = size (x_slope); if (m < n), x_slope = x_slope'; end
	[m,n] = size (y_slope); if (m < n), y_slope = y_slope'; end
	[m,n] = size (i_slope); if (m < n), i_slope = i_slope'; end
	[m,n] = size (j_slope); if (m < n), j_slope = j_slope'; end
	[m,n] = size (z_slope); if (m < n), z_slope = z_slope'; end
	n1 = length (z_slope);
end

% FRAM FRAM FRAM 
% scale and detrend 
if (n_args == 5)
    xsc = max(x_data)-min(x_data); 
    ysc = max(y_data)-min(y_data);
    x_data = x_data/xsc; x_out = x_out/xsc; 
    y_data = y_data/ysc; y_out = y_out/ysc; 
    [z_data,r0,gg] = detrend2(x_data,y_data,z_data);
elseif (n_args == 7)
    xsc = max([x_slope;x_out])-min([x_slope;x_out]);
    ysc = max([y_slope;y_out])-min([y_slope;y_out]);
    x_slope = x_slope/xsc; x_out = x_out/xsc;
    y_slope = y_slope/ysc; y_out = y_out/ysc;
    % scaling                  
    vx1=i_slope.*z_slope/xsc; 
    vy1=j_slope.*z_slope/ysc; 
    is1=vx1./sqrt(vx1.^2+vy1.^2);
    js1=vy1./sqrt(vx1.^2+vy1.^2);
    z_slope=sqrt(vx1.^2+vy1.^2);
    % detrending
    th1=acos(is1);
    phi=atan(gg(2)/gg(3));
    gam=th1-phi;    
    i_slope=cos(gam);
    j_slope=sin(gam);
else  % n_args == 10
    xsc = max(x_data)-min(x_data); 
    ysc = max(y_data)-min(y_data);
    x_data = x_data/xsc; x_out = x_out/xsc; x_slope=x_slope/xsc;
    y_data = y_data/ysc; y_out = y_out/ysc; y_slope=y_slope/ysc;    
    [z_data,r0,gg] = detrend2(x_data,y_data,z_data);
    % scaling 
    vx1=i_slope.*z_slope/xsc; %vx=i_slope.*z_slope; vx1=vx/xsc;
    vy1=j_slope.*z_slope/ysc; %vy=j_slope.*z_slope; vy1=vy/ysc;
    is1=vx1./sqrt(vx1.^2+vy1.^2);
    js1=vy1./sqrt(vx1.^2+vy1.^2);
    z_slope=sqrt(vx1.^2+vy1.^2);
    % detrending
    th1=acos(is1);
    phi=atan(gg(2)/gg(3));
    gam=th1-phi;    
    i_slope=cos(gam);
    j_slope=sin(gam);
end
% FRAM FRAM FRAM 
    
% Assembly final xp, yp, and zp vectors (possibly combination of data and slopes)

if (n_args == 10)	% z_data and z_slope supplied; put xyz in general point vector
	xp = [x_data; x_slope];
	yp = [y_data; y_slope];
	zp = [z_data; z_slope];
elseif (n_args == 7)	% z_slope supplied; put xyz in general point vector
	xp = x_slope;
	yp = y_slope;
	zp = z_slope;
else 			% z_data supplied; put xyz in general point vector
	xp = x_data;
	yp = y_data;
	zp = z_data;
end

% Now build the square n x n linear system that must be solved for the alpha's
if verbose
    disp ('build matrix')
    tic
end
n = n0 + n1;
A = zeros (n, n);
if (n_args == 5 | n_args == 10)	% z_data supplied; build data matrix rows
	for i = 1:n0
		r = (abs ((xp(i) - xp) + IM * (yp(i) - yp)));
		A(i,:) = (spline2dgreen (r, p))'; %FRAM, was spline2d_green
	end
end
if (n_args == 7 | n_args == 10)	% z_slope supplied; build slope matrix rows
	for i = 1:n1
		j = i + n0;
		dx = xp(j) - xp;
		dy = yp(j) - yp;
		A(j,:) = (spline2dgrad (dx, dy, i_slope(i), j_slope(i), p))'; %FRAM, was spline2d_grad 
	end
end
if verbose    
    disp ([num2str(toc),' now solve matrix'])
    tic
end
% Done building square linear system, now solve it

alpha = A \ zp;

if (nargout == 2)	% Return eigenvalues
	if verbose;disp ('find eigenvalues'); end
	l = svd (A);
end

if verbose
    disp ([num2str(toc),' now evaluate'])
    tic
end
% Now evaluate final solution at output locations

z = zeros (size(x_out));
for i = 1:length(alpha)
	r = abs ((x_out - xp(i)) + IM * (y_out - yp(i)));
    z = z + alpha(i) * (spline2dgreen (r, p)); %FRAM, was spline2d_green
end

if verbose    
    disp ([num2str(toc),' now done'])
end

% FRAM FRAM FRAM 
%retrend
z = z+gg(1)+(x_out-r0(1))*gg(2)+(y_out-r0(2))*gg(3);
