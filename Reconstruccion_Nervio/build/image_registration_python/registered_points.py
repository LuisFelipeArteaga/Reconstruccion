######################### Registro de Imagenes #################################

# Este script realiza el registro de las imagenes de profundidad, teniendo en cunta la posicion
# de cada imegen en una trayectoria alrededor de un objeto. La trayectoria se almacena en un archivo txt.
# Cada trayectoria  de posicion de cada imagen se almacena en un cuaternion (stamp,tx,ty,tz,qx,qy,qz,qw)
# Las imagenes se almacenan en formato png, el tamano de las imagenes es de 640x480 pixeles 

import numpy as np 
import argparse
import os
import sys
from functions import *


# Carpeta que contiene las imagenes de profundidad, para el registro de las imagenes.
# rgb.txt = almacena el nombre de las imagenes rgb
# depth.txt = almacena el nombre de la imagenes de profundidad
# trayectory.txt = almacena la trayectoria de la camara o sensor
dirc = "images/images/"
list_rgb = read_txt("images/images/rgb.txt")
list_depth = read_txt("images/images/depth.txt")
list_trajec = read_txt("images/trayectory/trayectory.txt")


# Como los archivos txt y las imegenes rgb y de profundidad no tienen el mismo nombre
# se empareja las imagenes con la lista de los txt se establece una diff minima(Los  nombre son
# numeros del tiempo que se tomaron)
diff1 = 0.02
diff2 = 0.01
# Empareja entre las imagenes de rgb y detph
matches1 = dict(match(list_rgb,list_depth,diff1))
# Empareja entre las imagnes de rgb y trayectoria
matches2 = match(matches1,list_trajec,diff2)
# matches = devuelve una lista que contiene el nombre de las imagenes y las  su trayectoria respectiva.
matches2.sort()


all_points = []

# De una lista original de imagenes se seleciona un numer menor de imagenes
number_images = int(sys.argv[1]) # Numero total de imagenes/nuber_images = Numero de imagenes para el registro


list = range(0,len(matches2),number_images)
# Recorre todas las imagenes de profundidad de las que se obtendran las nubes de puntos
for maxIter, i in enumerate(list):
	
	rgb_name, trajec_name = matches2[i]

	rgb_file = list_rgb[rgb_name][0]
	depth_file = list_depth[matches1[rgb_name]][0]
	pos = matrix_trans(list_trajec[trajec_name][0:8])

	# Se genera  una nube de puntos para cada imagen(Frame)
	downsample = int(sys.argv[2]) # submuestrea la nube de puntos.
	boundary = int(sys.argv[3]) # Se utiliza solo para tomar el objeto de interes de cada imagen de profundidad
	points = point_cloud(dirc+rgb_file,dirc+depth_file,pos,boundary,downsample)

	# Almacenar puntos de cada Frame
	all_points += points

# Ruta donde se almacena la nube de puntos
ruta = "point_cloud/"

# Almacena las nube de puntos en formato punto PLY.
out_filename = ruta+"noise_point_cloud.ply"
write_ply(out_filename,all_points)

