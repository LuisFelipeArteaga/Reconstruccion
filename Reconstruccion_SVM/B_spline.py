import numpy as np
import sys 
import os
import xml.etree.ElementTree as ET
from scipy import interpolate
from scipy import stats
import time

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Leer archivo xml, que contiene coordenadas x,y del contorno del nervio.
def read_xml(address):
	listXY = [];
	tree = ET.parse(address)
	root = tree.getroot()
	for pt in root.iter('pt'):
		x = pt.find('x').text
		y = pt.find('y').text
		x = int(x)
		y = int(y)
		listXY.append([x,y])
	_,ind = np.unique(listXY,axis = 0, return_index =True)

	return [listXY[ii] for ii in np.sort(ind)]
#Save txt 
def save_XYZ(point_cloud, name):
	file = open(name,'w')
	point_cloud = np.asarray(point_cloud)
	for xx in point_cloud:
		x,y,z = xx
		file.write('%s\t %s\t %s\t \n' %(x,y,z))
	file.close()

# Interpolacion de puntos (Lineal-Cubica)
def B_spline1(pixelsXY,points,n,degree):
	pixelsXY = np.array(pixelsXY)	
	ctr = np.concatenate((pixelsXY, [pixelsXY[0,:]]), axis=0) 
	xy = stats.zscore(ctr,axis = 0)
	x=xy[:,0] 
	y=xy[:,1]
	tck,u = interpolate.splprep([x,y],k= degree,s=0)
	u=np.linspace(0,1, num=points, endpoint=False)
	out = interpolate.splev(u,tck)
	z = n+np.zeros(out[0].shape[0])
	pXYZ = np.c_[z,out[0],out[1]]
	
	return pXYZ,x,y,out[0],out[1]


# Interpolacion de puntos (Lineal-Cubica)
def B_spline2(pixelsXY,points,n):
	degree = 3
	xy = stats.zscore(pixelsXY,axis = 0)
	xy = xy.tolist()
	xy = xy + xy[0:degree + 1]
	xy = np.array(xy)
	n_points = len(xy)
	x = xy[:,0]
	y = xy[:,1]
	t = range(len(x))
	ipl_t = np.linspace(1.0, len(xy) - degree, points)
	x_tup = interpolate.splrep(t, x, k=degree, per=1)
	y_tup = interpolate.splrep(t, y, k=degree, per=1)
	x_list = list(x_tup)
	xl = x.tolist()
	x_list[1] = [0.0] + xl + [0.0, 0.0, 0.0, 0.0]
	y_list = list(y_tup)
	yl = y.tolist()
	y_list[1] = [0.0] + yl + [0.0, 0.0, 0.0, 0.0]
	x_i = interpolate.splev(ipl_t, x_list)
	y_i = interpolate.splev(ipl_t, y_list)

	z = n+np.zeros(x_i.shape[0])
	pXYZ = np.c_[z,x_i,y_i]
	
	return pXYZ,x,y,x_i,y_i

#Plot inteporlacion de B_spline

def plot(x,y,x_i,y_i,Time):
	plt.figure(1)
	plt.clf()
	plt.plot(x, y , 'o-r', x_i, y_i, 'b*')
	plt.legend(['Points', 'Inter B-spline', 'True'],loc='best')
	plt.axis([min(x_i)-.5, max(x_i)+.5, min(y_i)-.5, max(y_i)+.5])
	plt.title('Slice '+ ii )
	
	plt.show()
	#plt.draw()
	#plt.pause(Time)

######################### Main #########################################
if __name__ == '__main__':
	
	'''address = sys.argv[1] ## Carpeta que contiene los contornos
	name_xyz = sys.argv[2] ## Nombre de la nube de puntos	
	n = float (sys.argv[3]) # Separacion entre slices
	grado = int (sys.argv[4]) # Grado de interpolacion
	points = int (sys.argv[5]) # Numero de puntos por contorno'''

	address1 = 'Ciatico/' # Carpeta que contiene las imagenes
	address2 = 'Ciatico/' # Carpeta que contiene las
	folders = os.listdir(('Imagenes/'+address1))
	for kk in folders:

		if kk == 'P02_ciatico01_Labeled':#kk.endswith('Labeled'):
			point_cloud = []
			name_folder = 'Imagenes/'+address1+kk
			name_xyz = 'Data/'+address2+kk+'_SEP.xyz'
			print name_folder

			Inter_L = True
			nn = .3 # Separacion entre slices
			degree = 1 # Grado de interpolacion
			points = 50 # Numero de puntos por contorno
			n = 0
			
			for _,_,files in os.walk(name_folder):
					files = [jj for jj in files if jj.endswith('.xml') ]
					print files
					files.sort(key = lambda item : int(item.split('.')[0]))
					frames = [int(jj.split('.')[0]) for jj in files]
					disF = [(frames[kk+1]-frames[kk])for kk in range(len(frames)-1)]
					disF.insert(0,0)
					print disF
					ll = 0
					for ii in files: 
						pixelsXY = read_xml((name_folder+'/'+ii))

						if 	disF[ll] > 20:
							n = nn*6 + n;
						else:
							n = nn*disF[ll] + n

						ll += 1
						if Inter_L:
							pXYZ,x,y,x_i,y_i = B_spline1(pixelsXY,points,n,degree)
							#plot(x,y,x_i,y_i,.5)
							point_cloud.extend(pXYZ)
						else:
							pXYZ,x,y,x_i,y_i = B_spline2(pixelsXY,points,n)	
							#plot(x,y,x_i,y_i,.5)
							point_cloud.extend(pXYZ)

			save_XYZ(point_cloud,name_xyz)   

			point_cloud = np.asarray(point_cloud)
			fig = plt.figure(2)
			ax = fig.add_subplot(111, projection='3d')
			ax.plot3D(point_cloud[:,0],point_cloud[:,1],point_cloud[:,2],'.r')
			plt.title(kk)
			plt.axis('equal')
			plt.draw()
			plt.pause(20)

