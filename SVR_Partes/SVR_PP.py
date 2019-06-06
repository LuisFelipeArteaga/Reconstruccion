#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 10:41:44 2019

@author: felipe
"""
from __future__ import division
import time
from joblib import dump, load

import numpy as np
import matplotlib.pyplot as plt 
from sklearn.svm import SVR
from mpl_toolkits.mplot3d import Axes3D 
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV

import nibabel as nib
import cv2
import open3d as opd
import os
#import pcl

from skimage import measure
import plotly.offline as py

from plotly.tools import FigureFactory as FF


############################################################################
def Feature_Scaling(X,a,b):
    return a + ((X-np.min(X))*(b-a))/(np.max(X)-np.min(X))

################################################################################
def load_point_cloud(sliceInf,sliceSup,path):
    #path = '/home/felipe/Desktop/Data_Heart/MM-WHS_2/mr_train/mr_train_03/Affine_ANTS_reg/mr_train_03Myc_seg.nii.gz'
    img  = nib.load(path)
    data = img.get_data()
    
    point_cloud = []
    if sliceSup == 'ALL':
        sliceSup = data.shape[0]
        sliceInf = 0
    
    for pp in np.linspace(sliceInf,sliceSup,data.shape[0],endpoint = False, dtype = int):
        sliceL = data[pp,:,:]
        #plt.figure('1')
        
        if 1 in sliceL:
            _, thresh = cv2.threshold(np.uint8(sliceL*255),127,255,cv2.THRESH_BINARY)
            contours,_ = cv2.findContours(thresh,cv2.RETR_TREE,cv2.CHAIN_APPROX_NONE)
            indA = []
            for jj in range(len(contours)):
                indA.append(cv2.contourArea(contours[jj]))
            ind = [ii for ii,jj in enumerate(indA) if jj<10]
            ind.sort(reverse=True)
            for ii  in ind:
                contours.pop(ii)
            if len(contours) != 0:
                indC = np.array(np.concatenate(contours))
                indC = np.reshape(indC,(len(indC),2))
                zpoints =  pp*np.ones((indC.shape[0],1))
                points = np.hstack((indC,zpoints))
                point_cloud.append(points)    

    points_cloud = np.concatenate(point_cloud,axis = 0)
    return points_cloud     
    #print(img.header)
######################Save PLY #########################
def save_points_cloud(ply_file,vertex):
    file = open(ply_file,"w")
    file.write('''ply
format ascii 1.0
element vertex %d
property float x
property float y
property float z
end_header
'''%len(vertex))
    for row in vertex:
        file.write('%f %f %f \n' %(row[0],row[1],row[2]))
  
    file.close()
    
#################################################################################
def save_points(name):
    vertex = load_point_cloud()
    vertexN = Feature_Scaling(vertex,-1,1)
    filename = os.path.join('Points_Clouds',name)
    save_points_cloud(filename,vertexN)
##################################################################################
def get_data_myc(name):
    
    filename = os.path.join('Points_Clouds',name)
    pcd  = opd.read_point_cloud(filename)
    
    print('compute normals ...')
    opd.estimate_normals(pcd, search_param = opd.KDTreeSearchParamHybrid(
            radius = .5, max_nn = 30))
    print('end')
    
    mycCenter = np.mean(np.asarray(pcd.points[:]),axis=0)
    for ii,qq in enumerate(pcd.points):
       p1 = mycCenter - qq;
       p2 = np.asarray(pcd.normals)[ii,:]
       angle = np.arctan2(np.linalg.norm(np.cross(p1,p2)),np.dot(p1,p2))
       if angle < np.pi/2 or angle < -np.pi/2:
           pcd.normals[ii] = -pcd.normals[ii]
    
    delta = 0.01;
    pointsOut = np.asarray(pcd.points) + delta*np.asarray(pcd.normals)
    pointsIn = np.asarray(pcd.points) - delta*np.asarray(pcd.normals)
    
    data = np.concatenate((pointsOut,np.asarray(pcd.points),pointsIn),axis = 0)
    labels = np.hstack((np.ones(pointsIn.shape[0],),
                        np.zeros(pointsIn.shape[0],),
                        -np.ones(pointsIn.shape[0],)))
    opd.draw_geometries([pcd])
    return data,labels

if __name__ == '__main__':
    for ii in range(9,14):

        path = '/home/felipe/Desktop/Data_Heart/MM-WHS_2/mr_train/mr_train_'+str(ii+1)+'/Affine_ANTS_reg/mr_train_'+str(ii+1)+'Myc_seg.nii.gz'
        sliceInf = 75
        sliceSup = 'ALL'
        mycPoints= load_point_cloud(sliceInf,sliceSup,path)
        mycPoints = Feature_Scaling(mycPoints,-1,1)
        save_points_cloud('Point_Cloud2/subj'+str(ii+1)+'.ply',mycPoints)
    