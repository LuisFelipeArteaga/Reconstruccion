function G = spline2dgreen (x, c)
%	$Id: spline2dgreen.m,v 1.4 2014/06/29 03:31:27 myself Exp $
% SPLINE2DT_GREEN	Green function for 2-D spline in tension
%
%	g = spline2dt_green (x, c)
%
%	x	- abscissa values
%	c	- tension parameter (c = sqrt (t/(1-t))
%
%  See Wessel, P, D. Bercovici, 1998, Gridding with Splines in Tension : A Green function Approach,
%       Math. Geol., 30, 77-93.

% Translated to Matlab from C to avoid needing to build mex executables.
% spline2d_green computes the Green function for a 2-d spline possibly
% in tension, G(u) = G(u) - log(u), where u = c * x and c = sqrt (t/(1-t)).
% The modified Bessel function K of order zero is based on Num. Rec.
% All x must be >= 0.  When c = 0 it degenerates to x^2 * log(x)
% 

if c == 0 % Just regular spline
	mask = find(x == 0);
	if length(mask)>0, x(mask) = exp(1)*ones(length(mask),1); end        
	G = (x.^2).*(log(x)-1.0);  % SaGa, proportional to Sandwell article.  
else    % In tension    
	ic = 1/c;
	g0 = 0.115931515658412420677337;	% log(2) - 0.5772156...
	mask = find(x == 0);
	if length(mask)>0, x(mask) = ones(length(mask),1); end        
	id = find(x<=2*ic);
	cx= c*x(id);
	t = cx.*cx;
	y = 0.25*t;
	z = t/14.0625;
	G(id) = (-log(0.5*cx) .* (z .* (3.5156229 + z .* (3.0899424 + z .* (1.2067492 + z .* (0.2659732 + z .* (0.360768e-1 + z .* 0.45813e-2))))))) + (y .* (0.42278420 + y .* (0.23069756 + y .* (0.3488590e-1 + y .* (0.262698e-2 + y .* (0.10750e-3 + y .* 0.74e-5))))));
	id = find(x>2*ic);    
	y = 2*ic./x(id);
	cx = c*x(id);
	G(id) = (exp (-cx) ./ sqrt (cx)) .* (1.25331414 + y .* (-0.7832358e-1 + y .* (0.2189568e-1 + y .* (-0.1062446e-1 + y .* (0.587872e-2 + y .* (-0.251540e-2 + y .* 0.53208e-3)))))) + log (cx) - g0;
	if length(mask)>0, G(mask) = zeros(length(mask),1); end        
	G=G';
end
