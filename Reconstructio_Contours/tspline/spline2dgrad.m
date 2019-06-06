function dgdn = spline2dgrad (dx, dy, ni, nj, c)
%	$Id: spline2dgrad.m,v 1.3 2014/06/29 03:31:27 myself Exp $
% SPLINE2DT_GRAD	Gradient of Green function for 2-D spline in tension
%
%	g = spline2dt_grad (x, c)
%
%	x	- abscissa values
%	c	- tension parameter (c = sqrt (t/(1-t))
%
%  See Wessel, P, D. Bercovici, 1998, Gridding with Splines in Tension : A Green function Approach,
%       Math. Geol., 30, 77-93.

% Translated to Matlab from C to avoid needing to build mex executables.
% spline2dgrad computes the radial derivative of Green function for a 2-d spline possibly
% in tension, G'(u) = 1/u - K1(u), where u = c * x and c = sqrt (t/(1-t)).
% The modified Bessel function K of order one is based on Num. Rec.
% All x must be >= 0.  When c = 0 it degenerates to x (2*log(x)-1)
%

if c == 0 	% Just regular spline
	%	r = hypot (dx(i), dy(i));
	r=sqrt(dx.^2+dy.^2);
	mask = find(r == 0);        
	if length(mask)>0,   r(mask)= ones(length(mask),1);end        
	dgdr = 2*log (r) - 1;  %Same as Sandwell [1984]
	if length(mask)>0,dgdr(mask)=zeros(length(mask),1);end
	dgdn = dgdr.*(dx.*ni+dy.*nj);
else    % In tension
	ic = 1/c;    
	r=sqrt(dx.^2+dy.^2);
	mask = find(r == 0);        
	if length(mask)>0, r(mask)= ones(length(mask),1);end        

	id=find(r <= 2*ic) ;
	cx = c*r(id);
	cx2 = cx.*cx;
	y = 0.25*cx2;
	z = cx2 / 14.0625;
	dgdr(id) = -((log(0.5.*cx).* (cx.* (0.5 + z.* (0.87890594 + z.* (0.51498869 + z.* (0.15084934 + z.* (0.2658733e-1 + z.* (0.301532e-2 + z.* 0.32411e-3)))))))) + (1.0./cx).* (y.* (0.15443144 + y.* (-0.67278579 + y.* (-0.18156897 + y.* (-0.1919402e-1 + y.* (-0.110404e-2 + y.* (-0.4686e-4)))))))); 
	dgdn(id) = dgdr(id) * (dx(id) * ni + dy(id) * nj);
	id=find(r > 2.0*ic) ;
	cx = c * r(id);
	icx = 1./cx;
	y = 2*icx;
	dgdr(id) = icx - ((exp (-cx)./ sqrt (cx)).* (1.25331414 + y.* (0.23498619 + y.* (-0.3655620e-1 + y.* (0.1504268e-1 + y.* (-0.780353e-2 + y.*(0.325614e-2 + y.* (-0.68245e-3))))))));
	dgdn(id) = c*dgdr(id).*((dx(id).* ni + dy(id).* nj)./r(id))';
	if length(mask)>0,dgdn(mask)=zeros(length(mask),1);end
	dgdn=dgdn';
end 
