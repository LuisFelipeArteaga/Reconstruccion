#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:24:34 2019

@author: felipe
"""

import numpy as np
import open3d as opd
import os
from numpy import linalg as LA
from sklearn.model_selection import KFold
from sklearn.svm import SVR
from joblib import dump
from joblib import load


##############################################################################
def crossval_val(X,t,nPart):
    kf = KFold(n_splits=nPart,shuffle=True,random_state=123)
    kf.get_n_splits(X)
    Train_Index = []
    Test_Index = []
    for train_index, test_index in kf.split(t):
        Train_Index.append(train_index)
        Test_Index.append(test_index)
    return Train_Index,Test_Index


##############################################################################
def crossval_val2(X,nPart):
    kf = KFold(n_splits=nPart,shuffle=True,random_state=123)
    kf.get_n_splits(X)
    Train_Index = []
    Test_Index = []
    for train_index, test_index in kf.split(X):
        Train_Index.append(train_index)
        Test_Index.append(test_index)
    return Train_Index,Test_Index

##############################################################################
def svrTrain(X,t,XTest,Kernel,sig,c,e):
     svr = SVR(kernel = 'rbf',gamma=sig, C=c, epsilon = e)
     svr.fit(X,t)
     svr.predict(XTest)
     return svr, svr.predict(XTest)
 
   
###############################################################################
def Feature_Scaling(X,a,b):
    return a + ((X-np.min(X))*(b-a))/(np.max(X)-np.min(X))
##############################################################################
def get_data_myc(filename):
 
    #filename = os.path.join('Points_Clouds',name)
    pcd  = opd.read_point_cloud(filename)
   
    auxPointsS = Feature_Scaling(np.asarray(pcd.points),-1,1)
    for pp,ii in enumerate(auxPointsS):
        pcd.points[pp] = ii
     
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
 
    delta = 0.0125;
    pointsOut = np.asarray(pcd.points) + delta*np.asarray(pcd.normals)
    pointsIn = np.asarray(pcd.points) - delta*np.asarray(pcd.normals)
 
    data = np.concatenate((pointsOut,np.asarray(pcd.points),pointsIn),axis = 0)
    labels = np.hstack((np.ones(pointsIn.shape[0],),
                        np.zeros(pointsIn.shape[0],),
                        -np.ones(pointsIn.shape[0],)))
    opd.draw_geometries([pcd])
    return data,labels   

###############################################################################
def Experimento_1(filename):
    return 0
###############################################################################

def Experimento_2(filename):
   
    Xdata,tLabels = get_data_myc(filename)
   
    Train_Index,Test_Index = crossval_val2(np.where(tLabels==1)[0],10)
    dimX = np.where(tLabels==1)[0].shape[0]
    ii = 0
    auxTrainIndex = Train_Index[ii]
    auxTestIndex = Test_Index[ii]
    Xtrain = Xdata[np.concatenate((auxTrainIndex,auxTrainIndex+dimX,auxTrainIndex+(2*dimX))),:]
    tTrain = tLabels[np.concatenate((auxTrainIndex,auxTrainIndex+dimX,auxTrainIndex+(2*dimX)))]
    Xtest = Xdata[np.concatenate((auxTestIndex,auxTestIndex+dimX,auxTestIndex+(2*dimX))),:]
    tTest = tLabels[np.concatenate((auxTestIndex,auxTestIndex+dimX,auxTestIndex+(2*dimX)))]
    return Xtrain, tTrain, Xtest, tTest

###############################################################################
if __name__ == '__main__':
     #path = '/home/felipe/Desktop/Data_Heart/MM-WHS_2/Validacion_Surface/Point_cloud'
    filename = os.path.join('Point_Cloud','Myc.ply')
    Xtrain, tTrain, Xtest, tTest = Experimento_2(filename)

    #--------------------------------------------------------
    #main loop
    Kernel ='rbf'
   
    roi = 'Superior'
    path1 = 'Experimento_2_' +roi
    c= np.load((path1+'/data_save/cDic.npy'))
    gamma = np.load((path1+'/data_save/gammaDic.npy'))
    epsilon = np.load((path1+'/data_save/epsilonDic.npy'))
    error = np.load((path1+'/data_save/errorDic.npy'))
    error2 = np.load((path1+'/data_save/error2Dic.npy'))
   
#   #--------------------------------------------
#   #Best Parametros
    #idxE2 = np.argmin(error2)
    for idxE2 in range(10):
        filenameSVR =os.path.join(path1,'Model','svrMyc'+str(idxE2+1)+'.joblib')
        print(filenameSVR)
        svrBest = load(filenameSVR)
       
        cbest = svrBest.C
        gammabest = svrBest.gamma
        epsilonbest = svrBest.epsilon
   
        print('C = ',cbest,' Gamma = ',gammabest,' Epsilon = ',epsilonbest)
        svr,t_est2 = svrTrain(Xtrain,tTrain,Xtest,Kernel,gammabest,cbest,epsilonbest)
        filename2 = os.path.join('Model_Experimento_2_'+roi,'svrMycSUP'+str(idxE2+1)+'.joblib')
        print('SAVE ',filename2,'\n')
        dump(svr, filename2)
   
   
        ## Error RMSD
        Error2 = np.zeros((1,2))
        Error2[0,0] =  np.sqrt((np.sum((t_est2-tTest)**2)/tTest.shape[0]))
        Error2[0,1] = LA.norm((t_est2 - tTest))/LA.norm(tTest)
        np.save('Model_Experimento_2_'+roi+'/ErrorMyc'+str(idxE2+1)+'.npy',Error2)