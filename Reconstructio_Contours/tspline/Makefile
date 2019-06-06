#	$Id: Makefile,v 1.1 2013/01/07 20:04:18 myself Exp $
# Makefile for tspline mex components.
# You will have to mess around to get the settings you need for your platform
#
MATDIR		= /usr/local/matlab4
CC		= /usr/local/SUNWspro/SC4.0/bin/cc
ARCH		= `arch`
EXT		= mexsol	# Or mexhp3 (HP/9000 s300), mexds (Dec), mexsg (SGI), mexvx (Vax), mexsol (Solaris)
#
#	The rest should should be ok
#
#MEX	= -L$(MATDIR)/extern/lib/$(ARCH) -lmex
MEX	= -L$(MATDIR)/extern/lib/sol2 -lmat 	# On Solaris
FLAGS	= -g -I.. -I$(MATDIR)/extern/include -v
CFLAGS	= -O -v

spline2d_green:	
		cmex $(FLAGS) spline2d_green.c $(MEX) 

spline2d_grad:	
		cmex $(FLAGS) spline2d_grad.c $(MEX) 

.c.o:
		$(CC) -c $(CFLAGS) $<
