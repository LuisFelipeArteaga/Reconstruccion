import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy.matlib as npm
import os


# Load archivo que contiene la nubes de puntos
def load_xyz(name):
    pointsC = np.loadtxt(name)

    return pointsC

# Rota un contorno en x.y and z
def rotation_cont(contour,rad):
    m = np.mean(contour[:,0])
    contour[:,0] = contour[:,0]- m
    theta1 = -rad +(rad+rad)*np.random.rand(1)    
    theta2 = -rad +(rad+rad)*np.random.rand(1)
    theta3 = -rad +(rad+rad)*np.random.rand(1)
    Rx = [[1,0,0],[0,np.cos(theta1),-np.sin(theta1)],[0,np.sin(theta1),np.cos(theta1)]]
    Ry = [[np.cos(theta2),0,np.sin(theta2)],[0,1,0],[-np.sin(theta2),0,np.cos(theta2)]]
    Rw = [[np.cos(theta3),-np.sin(theta3),0],[np.sin(theta3),np.cos(theta3),0],[0,0,1]]

    R1 = np.dot(Rx,contour.T) 
    R2 = np.dot(Ry,R1)
    R_contour = np.dot(R2.T,Rw)
    R_contour[:,0] = R_contour[:,0] + m
    
    return R_contour

 # Calcula la normal al plano
def normal_plane(R_contour):
    posP = np.round(R_contour.shape[0]/3)
    P = R_contour[1,:]
    Q = R_contour[posP,:]
    R = R_contour[posP*2,:]
    PQ = Q-R
    PR = R-P
    n = np.cross(PQ,PR)
    n_plane = n/np.linalg.norm(n)

    return n_plane


def normal_contour(R_contour,n_plane):
    dpx = np.gradient(R_contour[:,0])
    dpy = np.gradient(R_contour[:,1])
    dpz = np.gradient(R_contour[:,2])
    n_contour = np.cross(np.c_[dpx,dpy,dpz],n_plane)
    vecU = npm.repmat(np.linalg.norm(n_contour,axis = 1),3,1)
    n_contour = np.divide(n_contour,vecU.T)
    sensorCenter = np.mean(R_contour)

    #p1 = sensorCenter - R_contour
    return n_contour

def plot_quiver(x,y,z,u,v,w,Time,jj):
    fig = plt.figure(1)
    plt.clf()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot3D(x,y,z,'.r')
    ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True)
    plt.title('Normals contour '+jj)
    ax.view_init(elev=180., azim= 270)
    plt.show()
    #plt.draw()
    #plt.pause(Time)

def plot_quiver1(x,y,u,v,Time,jj):
    plt.figure(1)
    plt.clf()
    plt.plot(x,y,'.r')
    plt.quiver(x, y, u, v,  units='x', pivot='tail', width=0.022, scale=1 / 0.15)
    #plt.title('Normals contour '+jj)
    plt.axis([min(x)-.5, max(x)+.5, min(y)-.5, max(y)+.5])
    plt.axis('off')
    plt.show()



    
    #plt.pause(Time)

def labels_points(point_cloud, normals,delta):
    label1 = np.zeros((len(point_cloud),1))
    label2 = np.ones((len(point_cloud),1))
    point_cloud = np.asarray(point_cloud)

    np.savetxt('normals.txt', normals, delimiter=' ') 
    np.savetxt('point_cloud.txt', point_cloud, delimiter=' ') 

    sur = np.c_[label1,point_cloud]
    outside = np.c_[label2,point_cloud+(delta*np.array(normals))]
    inside = np.c_[-label2,point_cloud-(delta*np.array(normals))]
    
    return np.r_[sur,outside,inside] 

def save_point_Labels(point_cloud_L, name):
    file = open(name,'w')
    point_cloud_L = np.asarray(point_cloud_L)
    for xx in point_cloud_L:
        L,x,y,z = xx
        file.write('%s 1:%s\t 2:%s\t 3:%s\t \n' %(L,x,y,z))
    file.close()   

######################### Main #########################################

if __name__ == '__main__':

        
    name = 'Data/data_HVSMR.xyz'
    pointsC = load_xyz(name)
    _,indices = np.unique(pointsC[:,0], return_index=True)
    if pointsC.shape[1] == 4:
        pointsC = pointsC[:,1:4]
    print pointsC.shape    
    point_cloud = []
    normals = []
    rot = 0
    for ii in range(len(indices)-1):
        R_contour = pointsC[indices[ii]:indices[ii+1],:]
        #R_contour = rotation_cont(contour,rot)
        n_plane = normal_plane(R_contour) 
        n_contour = normal_contour(R_contour,n_plane)       

        point_cloud.extend(R_contour)
        normals.extend(n_contour)
        
        x = R_contour[:,1]
        y = R_contour[:,2]
        u = n_contour[:,1]
        v = n_contour[:,2]
        #plot_quiver1(x,y,u,v,.5,'d')

        

    name_data = 'Normals/data_HVSMR.xyz'
    if rot == 0: 
        print name_data

        save_point_Labels(labels_points(point_cloud, normals,0.05),(name_data+'.data'))
    else :
        save_point_Labels(labels_points(point_cloud, normals,0.05),(name_data+'_10G.data'))
    point_cloud = np.asarray(point_cloud)
    normals = np.asarray(normals)

    x = point_cloud[:,0]
    y = point_cloud[:,1]
    z = point_cloud[:,2]
    u = normals[:,0]
    v = normals[:,1]
    w = normals[:,2]
    #plot_quiver(x,y,z,u,v,w,.5,'d')  
    
    