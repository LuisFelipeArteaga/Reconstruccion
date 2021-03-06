#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 16:22:06 2019

@author: felipe
"""
##############
#Generar Surface

import numpy as np
import time 
import matplotlib.pyplot as plt 
import os

from skimage import measure
import plotly.offline as py

from joblib import load

from plotly.tools import FigureFactory as FF

######################Save PLY #########################
def write_ply(ply_file,vertex,faces):
  file = open(ply_file,"w")
  file.write('''ply
format ascii 1.0
element vertex %d
property float x
property float y
property float z
element face %d
property list uchar int vertex_index
end_header
'''%(len(vertex),len(faces)))

  for row in vertex:
    file.write('%f %f %f \n' %(row[0],row[1],row[2]))
  for row in faces:
    file.write('3 %i %i %i \n' %(row[0],row[1],row[2]))
  file.close()
  
###############################################################################
def marching_cubes(surf_eq, aspR,isovalor,name,name_folder):
    verts, faces, normals, values = measure.marching_cubes_lewiner(surf_eq, isovalor) 
    x,y,z = zip(*verts)
    print(type(verts))
    y1 = np.array(y)*aspR
    y = tuple(y1)
    z1 = np.array(z)*aspR
    z = tuple(z1)
    
    fig = FF.create_trisurf(x=x,
                            y=y, 
                            z=z, 
                            plot_edges=True,
                            show_colorbar=False,
                            showbackground=False,gridcolor='rgb(0,0,0)',
                            simplices=faces,
                            title="Isosurface",
                            aspectratio=dict(x=1, y=1, z=1))
    #print len(verts) ,len(faces)
    print(verts[:,0].shape)
    verts2 = np.c_[verts[:,0],y1,z1]
    name_folder1 = 'Surface_'+name_folder+'/' 
    write_ply(name_folder1+name+'.ply',verts2,faces)
    py.plot(fig,filename = ('Surface_Myc/'+name))
    
###############################################################################
def Feature_Scaling(X,a,b):
    return a + ((X-np.min(X))*(b-a))/(np.max(X)-np.min(X))

###############################################################################
def regresion(svr,n):

    print('implicit Surface ...')

    delta = 0.01
    xmin = -1 - delta
    xmax = 1 + delta
    ymin = -1 - delta
    ymax = 1 + delta
    zmin = -1 - delta
    zmax = 1 + delta
    
    X,Y,Z = np.mgrid[xmin:xmax:n*1j, ymin:ymax:n*1j, zmin:zmax:n*1j]
    XX = np.vstack((X.ravel(),Y.ravel(),Z.ravel())).T
    print(XX.shape)
    IS = svr.predict(XX)
    print(IS.shape)
    return IS

###############################################################################
if __name__ == '__main__':
    #name_folder = 'Experimento_1_Inferior'
    #name_model = 'svrMycINF2'
    
    name_folder = 'Experimento_2_Superior'
    slices  = 50
        
    for ii in range(1):
        name_model = 'svrMycSUP' + str(ii+1)
        filenameModel = os.path.join('Model_'+name_folder,name_model+'.joblib')
        svr = load(filenameModel)
        #-----------------------------------------------------------------------
        t0 = time.time()
        print('compute regresion ....')
        f = regresion(svr,slices)
        print(np.min(f))
        print(np.max(f))
        f = f.reshape([slices,slices,slices])
        reg_fit = time.time() - t0
        print('fitted in %.3f s \n\n' % reg_fit)

        for ii in np.linspace(0,slices-1,7,dtype= int):
            plt.figure('Slice regresion')
            plt.imshow(f[:,:,ii])
            plt.show()
        
        #------------------------------------------------------------------------
        isovalor = 0.1
        aspR = 1.0
        name2 = name_model
        t0 = time.time()
        print('compure surface ...')
        marching_cubes(f,aspR,isovalor,name2,name_folder)
        sur_fit = time.time() - t0
        print('fitted in %.3f s \n\n' % sur_fit)