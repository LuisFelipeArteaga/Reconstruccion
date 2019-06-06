import plotly.offline as py
import plotly.graph_objs as go
from plotly.tools import FigureFactory as FF
#from mayavi import mlab 
from plyfile import PlyData, PlyElement

#py.tools.set_credentials_file(username='lufarteagada', api_key='3if6FRpEaMT6GzEQvfnx')

import numpy as np
from skimage import measure
import os
import sys


############## Train svm ############################

def train_svm(name_data,options,name_model):

    print 'Compute SVM ...'
    comand2 = 'libsvm-3.22/svm-train'+' '+options+' '+name_data+' '+name_model
    print comand2
    os.system(comand2)

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



################ Marching Cubes ######################

def marching_cubes(surf_eq, aspR,isovalor,name,folder):
    verts, faces, normals, values = measure.marching_cubes_lewiner(surf_eq, isovalor) 
    x,y,z = zip(*verts)
    print type(verts)
    y1 = np.array(y)*aspR
    y = tuple(y1)
    z1 = np.array(z)*aspR
    z = tuple(z1)
    

    colormap=['rgb(0,0,255)','rgb(0,255,0)']


    fig = FF.create_trisurf(x=x,
                            y=y, 
                            z=z, 
                            plot_edges=False,
                            show_colorbar=False,
                            showbackground=False,gridcolor='rgb(0,0,0)',
                            simplices=faces,
                            title="Isosurface",
                            aspectratio=dict(x=1, y=1*aspR, z=1*aspR))
    #print len(verts) ,len(faces)
    print verts[:,0].shape
    verts2 = np.c_[verts[:,0],y1,z1]
    write_ply('Surface/'+folder+name+'.ply',verts2,faces)
    py.plot(fig,filename = ('Surface/'+folder+name))
 

############### Parametros svm ###########################

def get_svm(name_model):
    i = 0
    w = []
    xi = []
    print 'Get_svm ...'
    Num_feat = 0
    for line in open(name_model):
        line  = line.split(None, 1)
        if len(line) == 1: line += ['']
        label , features = line
        if label == 'gamma':
            gamma = float(line[1])
        if label == 'total_sv':
            total_sv = int(line[1])
        if label == 'rho':
            rho = float(line[1])
        if i == 1:
            val = features.replace(':',' ').split(' ')
            w += [float(label)]
            xi += [np.array(val[1::2],dtype = float)]
            if Num_feat == 0:
                Num_feat = int(val[-3])
        if label == 'SV':   
            i = 1

    w = np.asarray(w, dtype = float)
    xi = np.asarray(xi) 
    return Num_feat,gamma,total_sv,rho,w,xi

############ Compure Iso-surface ##################

def implicit_suface (name_model,n):

    Num_feat,gamma,total_sv,rho,w,xi = get_svm(name_model)

    print 'implicit Surface ...'

    xmin = np.min(xi[:,0]) - .25
    xmax = np.max(xi[:,0]) + .25
    ymin = np.min(xi[:,1]) - .25
    ymax = np.max(xi[:,1]) + .25
    zmin = np.min(xi[:,2]) - .25
    zmax = np.max(xi[:,2]) + .25

    xdis = np.absolute(xmin-xmax)
    ydis = np.absolute(ymin-ymax)
    zdis = np.absolute(zmin-zmax)

    
    X,Y,Z = np.mgrid[xmin:xmax:n*1j, ymin:ymax:n*1j, zmin:zmax:n*1j]
    auxOnes = np.ones((n,n,n))
    f = np.zeros((n,n,n))
    norm = np.zeros((n,n,n))

    for ii in range(w.shape[0]):
        norm = (X-xi[ii,0]*auxOnes)**2 + (Y-xi[ii,1]*auxOnes)**2 + (Z-xi[ii,2]*auxOnes)**2 
        f += w[ii]*np.exp(-gamma*np.absolute(norm))
        print w[ii]
    f = f-rho
    return f , xdis ,ydis, zdis

'''def implicit_suface (name_model,n):

    Num_feat,gamma,total_sv,rho,w,xi = get_svm(name_model)

    print 'implicit Surface ...'

    xmin = np.min(xi[:,0]) - .25
    xmax = np.max(xi[:,0]) + .25
    ymin = np.min(xi[:,1]) - .25
    ymax = np.max(xi[:,1]) + .25
    zmin = np.min(xi[:,2]) - .25
    zmax = np.max(xi[:,2]) + .25

    xdis = np.absolute(xmin-xmax)
    ydis = np.absolute(ymin-ymax)
    zdis = np.absolute(zmin-zmax)

    
    X,Y,Z = np.mgrid[xmin:xmax:n*1j, ymin:ymax:n*1j, zmin:zmax:n*1j]
    auxOnes = np.ones((n,n,n))
    f = np.zeros((n,n,n))
    norm = np.zeros((n,n,n))

    for ii in range(w.shape[0]):
        prod = (X*xi[ii,0]) + (Y*xi[ii,1]) + (Z*xi[ii,2]) 
        f += w[ii]*((gamma*prod + 10)**3)
    f = f-rho

    return f , xdis ,ydis, zdis'''

##################### Main #########################


if __name__ == '__main__':

    folder = 'Ciatico/'

    files_Labels = os.listdir('Normals/'+folder)

    for ii in files_Labels:  
        if ii == 'P02_ciatico01_Labeled_SEP.data':#ii.endswith('10G.data'):    
            name =ii[0:-5]
            print '\n\n',name,'\n\n'
            name_model = 'Model/'+folder+name+'.model'

            gamma = 1/((.9)**2)


            print gamma
            param = '-s 3 -t 2 -c 1 -g '+str(gamma)+' -p 0.1' # Kernel Gaussiano exp(-gamma*|u-v|^2)
            #param = '-s 3 -t 1 -c 1 -r 10 -g .2 -d 3 ' #Kernel Polinomial (gamma*u'*v + coef0)^degree

            #train_svm('Normals/'+folder+name+'.data',param,name_model)

            isovalor = 0.0
            slices  = 100
            f , xdis ,ydis, zdis = implicit_suface(name_model,slices)
            aspR = np.mean([ydis,zdis])/ xdis
            name2 = 'surface4'

            marching_cubes(f,aspR,isovalor,name2,folder)