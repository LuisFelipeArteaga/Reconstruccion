#!/usr/bin/python
# Software License Agreement (BSD License)
#
# Copyright (c) 2013, Juergen Sturm, TUM
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above
#    copyright notice, this list of conditions and the following
#    disclaimer in the documentation and/or other materials provided
#    with the distribution.
#  * Neither the name of TUM nor the names of its
#    contributors may be used to endorse or promote products derived
#    from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.





import numpy as np 
import argparse
import os
import sys
from PIL import Image
import struct 

# Parametros de la camara.
focalLength = 525.0
centerX = 319.5
centerY = 239.5
scalingFactor = 5000.0
_EPS = np.finfo(float).eps * 4.0

# Lee los nombres de las imagenes RGB, DEPTH y Tracyectorias
# Almacena los datos en un diccionario 
def read_txt(filename):
	name = open(filename)
	data = name.read()
	change = data.replace(","," ").replace("\t"," ").split("\n") 
	list = [[v.strip() for v in line.split(" ") if v.strip()!=""] for line in change if len(line)>0 and line[0]!="#"]
	list = [(float(l[0]),l[1:]) for l in list if len(l)>1]
	return dict(list)

# Empajera los nombres de las imagenes RGB y Depth, asi como el nombre de las imagenes con 
# la posicion  trayectoria definida
def match(list1, list2, difference):
	list1_key = list1.keys();
	list2_key = list2.keys();

	matches_p = [(abs(p - q), p ,q) for p in list1_key for q in list2_key if abs(p-q) < difference] 
	#print matches_p
	matches_p.sort();
	matches = [];
	for resta, p, q, in matches_p:
		if p in list1_key and q in list2_key:
			list1_key.remove(p)
			list2_key.remove(q)
			matches.append((p,q))
	matches.sort()
	return matches

# Transforma el vector de posciones (stamp,tx,ty,tz,qx,qy,qz,qw) 
# en una matriz de transformacion rigida
def matrix_trans(v):
	pos = v[0:3]
	quat = np.array(v[3:7],dtype=np.float64)
	dquat = np.dot(quat,quat);
	if dquat < _EPS:
		m1 = np.array(((1.0,0.0,0.0,pos[0])(0.0,1.0,0.0,pos[1])
			        (0.0,0.0,1.0,pos[2])(0.0,0.0,0.0,1.0)), dtype=np.float64)
		return m1
	quat *= np.sqrt(2.0/dquat)
	quat = np.outer (quat,quat)
	m2 = np.array(((1.0-quat[1,1]-quat[2,2],quat[0,1]-quat[2,3],quat[0,2]+quat[1,3],pos[0]),
        		(quat[0,1]+quat[2, 3],1.0-quat[0,0]-quat[2,2],quat[1,2]-quat[0,3], pos[1]),
        		(quat[0, 2]-quat[1, 3],quat[1, 2]+quat[0, 3], 1.0-quat[0, 0]-quat[1, 1],pos[2]),
        		( 0.0,0.0,0.0, 1.0)),dtype=np.float64)
	return m2

#Genera la nube de puntos para cada imagen a partir de la imagen de profundidad y rgb

def point_cloud (rgb, depth, pos,limit, downsample):
	img_rgb = Image.open(rgb)
	img_depth = Image.open(depth)
	# Recorre toda la imegen de profundidad
	points = []
	for v in range(0,img_rgb.size[1],downsample):
	    for u in range(0,img_rgb.size[0],downsample):
	        # Se extraen los pixeles de cada imagen
	        color = img_rgb.getpixel((u,v))
	        Z = img_depth.getpixel((u,v))
	        # Se extablece un limite para tomar solo el objeto de interes
	        if Z != 0 and Z < limit:
	        	# Se genera una punto 3D para cada pixel y se ubica con respecto a la posicion de cada imagen
	        	Z = Z/scalingFactor
	        	X = (u-centerX)*Z/focalLength
	        	Y = (v-centerY)*Z/focalLength
	        	vec = np.matrix([[X],[Y],[Z],[1]])
	        	vec_transf = np.dot(pos,vec) # Coordenadas del punto.
		        points.append("%f %f %f %d %d %d 0\n"%(vec_transf[0,0],vec_transf[1,0],vec_transf[2,0],0.0,88.0,241.0))
	
	return points

#Almacena la nube de puntos en formato punto PLY
def write_ply(ply_file,points):
    file = open(ply_file,"w")
    file.write('''ply
format ascii 1.0
element vertex %d
property float x
property float y
property float z
property uchar red
property uchar green
property uchar blue
property uchar alpha
end_header
%s
'''%(len(points),"".join(points)))
    file.close()
    #print "Nube de puntos almaceda %s tamano %f"%(ply_file,len(points))

