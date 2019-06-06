import numpy as np 
import sys 
import os
import xml.etree.ElementTree as ET
import cv2
from scipy.ndimage.filters import gaussian_filter

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy import stats



print '\n\nGenera la nube de puntos de las imagenes segmentadas '
print 'python Get_point_cloud.py address  name.xyz SepZ plot3d'
print 'address = carpeta que contiene los archivos txt y xml'
print 'name.xyz = nombre del archivo donde se almacena la nube de puntos'
print 'SpeZ8(float) = Separacion entre slices del nervio'
print 'plot3d =  True o False \n\n'

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
	return listXY	

#Leer archivo txt, que contiene coordenadas x,y del contorno del nervio.
def read_txt(address):
	return np.loadtxt(address)

# Grafica la nube de puntos
def plot3d(point_cloud):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot3D(point_cloud[:,0],point_cloud[:,1],point_cloud[:,2],'.r')
	plt.show()

#Save txt 
def save_XYZ(point_cloud, name):
	file = open(name,'w')
	for xx in point_cloud:
		x,y,z = xx
		file.write('%s\t %s\t %s\t \n' %(x,y,z))
	file.close()
		
if __name__ == '__main__':
	
	# Ruta archivos txt o xml
	address = sys.argv[1]
	point_cloud = []
	for root,_,files in os.walk(address):
		files = [jj for jj in files if jj.endswith('.txt') | jj.endswith('.xml') ]
		files.sort(key = lambda item : int(item.split('.')[0]))
		
		n = 0;
		for ii in files: 
			# Matrix de ceros del size de la imagenes
			matrix = np.zeros((360,279,3), dtype = 'uint8')
			if ii.endswith('.txt'):
				pixelsXY = read_txt((address + ii))
				pixelsXY = [pixelsXY.astype(np.int32)]
			elif ii.endswith('.xml'):
				pixelsXY = read_xml((address + ii))	
				pixelsXY = np.array(pixelsXY)	

			# Extrae los contornos y centrar los datos 
			if len(pixelsXY) > 2: 
				cv2.drawContours(matrix,[pixelsXY],1,(255,0,0),-1)
				plt.imshow(matrix)
				plt.plot(pixelsXY[:,0],pixelsXY[:,1],'.b')
				plt.show()

				pX , pY = np.where(matrix[:,:,0] >= 250 )
				pX = pX - np.mean(pX)
				pY = pY - np.mean(pY)
				pZ = np.zeros(pX.shape[0])
				pXYZ = np.c_[n+pZ,pX,pY]

				n += float (sys.argv[3])
				point_cloud.extend(pXYZ)


	point_cloud = np.asarray(point_cloud)
	print point_cloud.shape
	# Dejar los datos entre -1 and 1
	point_cloud[:,1:4] = stats.zscore(point_cloud[:,1:4],axis = 0)
	#dataG = gaussian_filter(point_cloud[:,1:4] , 0.5)

	save_XYZ(point_cloud,sys.argv[2])
	if bool(sys.argv[4]):			
		plot3d(point_cloud)
	
			

		
