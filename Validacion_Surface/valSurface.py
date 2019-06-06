#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 09:33:33 2019

@author: felipe
"""

import numpy as np
#import matplotlib.pyplot as plt
#import time


#import nibabel as nib
#import cv2
import open3d as opd
import os

import itertools
from functools import partial
from multiprocessing import Pool

from sklearn.model_selection import KFold
#from scipy import stats
from scipy import spatial
from sklearn.svm import SVR
#from mpl_toolkits.mplot3d import Axes3D
#from sklearn.preprocessing import StandardScaler
#from sklearn.model_selection import GridSearchCV
from numpy import linalg as LA
from joblib import dump 
#from loblib import load

##############################################################################
def crossval_bal(X,t,nPart):
    kf = KFold(n_splits=nPart,shuffle=True,random_state=123)
    kf.get_n_splits(X)
    Train_Index = []
    Test_Index = []
    for train_index, test_index in kf.split(t):
        Train_Index.append(train_index)
        Test_Index.append(test_index)
    return Train_Index,Test_Index

###############################################################################
def Feature_Scaling(X,a,b):
    return a + ((X-np.min(X))*(b-a))/(np.max(X)-np.min(X))

##############################################################################
def svrTrain(X,t,XTest,Kernel,sig,c,e):
     svr = SVR(kernel = 'rbf',gamma=sig, C=c, epsilon = e)
     svr.fit(X,t)
     svr.predict(XTest)
     return svr, svr.predict(XTest)
 
##############################################################################
def svrTrainP(X,t,XTest,tTest,param):
     sig = param[0]
     c = param[1]
     e = param[2]
     svr = SVR(kernel = 'rbf',gamma=sig, C=c, epsilon = e)
     svr.fit(X,t)
     t_est = svr.predict(XTest)
     Error = LA.norm((t_est - tTest))/LA.norm(tTest)
     outPutParam = [Error,sig,c,e]
     return outPutParam
  
##############################################################################
def FindParamSVR(X,t,Kernel,t_data,ii):
  
    #
    num_Workers = 4
    nP = 10
    Train_Index2,Test_Index2 = crossval_bal(X,t,nP)
    #-------------------------------------------------
    #exp1 = [-9, -7, -5, -3, -1, 0, 1, 3, 5, 7, 9]
    exp1 = [-9, -3, -1, 0, 1, 3, 9]
    c = np.exp2(exp1)
    epsilon =  [0.0,0.0005,0.01,0.05,0.1,0.2]
  
    if Kernel == 'rbf':
        Nsig = 10
        s0 = np.median(spatial.distance.pdist(X,metric='euclidean'))
        sig = np.linspace(0.1*s0,s0,Nsig)
    else:
        sig = 1
        Nsig = 1
    #----------------------------------------------
    ## Partition
    cList = []
    gammaList = []
    epsilonList = []
    errorList = []

   
    for ff in range(nP):
        Xtrain = X[Train_Index2[ff]]
        Xtest = X[Test_Index2[ff]]
      
        tTrain = t[Train_Index2[ff]]
        tTest = t_data[Test_Index2[ff]]
      
        #------------------------------------------------------------
        param= list(itertools.product(sig,c,epsilon))
        #-------------------------------------------------------------
        with Pool(num_Workers) as p :
            funcSVR = partial(svrTrainP,Xtrain,tTrain,Xtest,tTest)
            Output = p.map(funcSVR,param)
          
        #---------------------------------------------------------
        Output = np.asarray(Output)
        minArg = np.argmin(Output[:,0])
        print(Output[minArg,1],Output[minArg,2],Output[minArg,3])
      
        gammaList.append(Output[minArg,1])
        cList.append(Output[minArg,2])
        epsilonList.append(Output[minArg,3])
        errorList.append(Output[minArg,0])
   
        filename2 = os.path.join('data_save','Error','error'+str(ii)+'_'+str(ff)+'.npy')
        np.save(filename2,cDic)
    return cList, gammaList, epsilonList, errorList
              
          
##############################################################################
def get_data_myc(filename):
  
    #filename = os.path.join('Points_Clouds',name)
    pcd1  = opd.read_point_cloud(filename)
    
    auxPointsS = Feature_Scaling(np.asarray(pcd1.points),-1,1)
    for pp,ii in enumerate(auxPointsS):
        pcd1.points[pp] = ii
      
    pcd = opd.voxel_down_sample(pcd1, voxel_size = 0.06)
    #pcd=pcd1
    print('Point Cloud size =  ',pcd.dimension)
  
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
      
###############################################################################
###############################################################################

if __name__ == '__main__':
    #path = '/home/felipe/Desktop/Data_Heart/MM-WHS_2/Validacion_Surface/Point_cloud'
    nPart = 10
    filename = os.path.join('Point_Cloud','bunny.ply')
    X,t = get_data_myc(filename)
   
    #---------------------------------------------------------
    #Validation
    Train_Index,Test_Index = crossval_bal(X,t,nPart)
    print(X.shape[0]/3)
    #---------------------------------------------------------
    #Classification initialization
    #--------------------------------------------------------
    #main loop
    Kernel ='rbf'
    cDic = {}
    gammaDic = {}
    epsilonDic = {}
    errorDic = {}
  
    errorList2 = []
      
    for ii in range(1):
        print('Fold ',ii+1,' of ' ,nPart)
        #-----------------------------------------------------
        #Partition
        Xtrain = X[Train_Index[ii]]
        Xtest = X[Test_Index[ii]]
        #Xtest = np.reshape(Xtest,(Xtest.shape[0],1))
      
        tTrain = t[Train_Index[ii]]
        tTest = t[Test_Index[ii]]
#      
        y_dataR = t[Train_Index[ii]]
#      
#
        cList, gammaList, epsilonList, errorList = FindParamSVR(
                                                    Xtrain,tTrain,Kernel,y_dataR,ii)
#        #--------------------------------------------
#        #Best Parametros
        indx = np.argmin(errorList)
#      
        cbest = cList[indx]
        gammabest = gammaList[indx]
        epsilonbest = epsilonList[indx]

        svr,t_est2 = svrTrain(Xtrain,tTrain,Xtest,Kernel,gammabest,cbest,epsilonbest)
        filename2 = os.path.join('Model','svrBunny'+str(ii+1)+'.joblib')
        print('SAVE ',filename2,'\n')
        dump(svr, filename2)
        Error2 = LA.norm((t_est2 - tTest))/LA.norm(tTest)
      
        errorList2.append(Error2)
#      
#        #-------------------------------------------
#        #save
        cDic[str(ii)] = cList
        gammaDic[str(ii)] = gammaList
        epsilonDic[str(ii)] = epsilonList
        errorDic[str(ii)] = errorList
#      
        np.save('data_save/cDicBunny.npy',cDic)
        np.save('data_save/gammaDicBunny.npy',gammaDic)
        np.save('data_save/epsilonDicBunny.npy',epsilonDic)
        np.save('data_save/errorDicBunny.npy',errorDic)
        np.save('data_save/error2DicBunny.npy',errorList2)