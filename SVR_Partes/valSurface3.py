#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 09:33:33 2019

@author: felipe
"""

import numpy as np
import open3d as opd
import os

import itertools
from functools import partial
from multiprocessing import Pool

from sklearn.model_selection import KFold
from scipy import spatial
from sklearn.svm import SVR
from joblib import dump 


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
     ## Error RMSD
     Error =  np.sqrt((np.sum((t_est-tTest)**2)/tTest.shape[0]))
     ##
     #Error = LA.norm((t_est - tTest))/LA.norm(tTest)
     outPutParam = [Error,sig,c,e]
     return outPutParam
  
##############################################################################
def FindParamSVR(X,t,Kernel,t_data,gg):
  
    #
    num_Workers = 4
    nP = 10
    Train_Index2 = []
    Test_Index2 = []
    tTrain_Index2 = []
    tTest_Index2  = []
    Train_IndexAux,Test_IndexAux = crossval_val2(t_data,nP)
    for ii in Train_IndexAux:
        #print(ii)
        auxPos = np.concatenate((ii,ii+t_data.shape[0],ii+(2*t_data.shape[0])))
        Test_Index2.append(ii+t_data.shape[0])
        Train_Index2.append(auxPos)
    for jj in Test_IndexAux:
        #print(jj)
        auxPos = np.concatenate((jj,jj+t_data.shape[0],jj+(2*t_data.shape[0])))
        tTest_Index2.append(jj+t_data.shape[0])
        tTrain_Index2.append(auxPos)
    #-------------------------------------------------
    exp1 = [-9, -7, -5, -3, -1, 0, 1, 3, 5, 7, 9].
    c = np.exp2(exp1)
    epsilon =  [0.001,0.005,0.015,0.025,0.1,0.2]
    #epsilon =  [0.0005,0.015,0.1]
 
    Nsig = 10
    s0 = np.median(spatial.distance.pdist(X,metric='euclidean'))
    sig = np.linspace(0.1*s0,s0,Nsig)

        
    #print(s0)
    #print(sig)
    #----------------------------------------------
    ## Partition
    cList = []
    gammaList = []
    epsilonList = []
    errorList = []

   
    for ff in range(nP):
        Xtrain = X[Train_Index2[ff]]
        Xtest = X[tTest_Index2[ff]]
      
        tTrain = t[Train_Index2[ff]]
        tTest = t[tTest_Index2[ff]]
        print('summa ',np.sum(tTest),Xtrain.shape,Xtest.shape,tTrain.shape,tTest.shape)
      
        #------------------------------------------------------------
        param= list(itertools.product(sig,c,epsilon))
        #-------------------------------------------------------------
        with Pool(num_Workers) as p :
            funcSVR = partial(svrTrainP,Xtrain,tTrain,Xtest,tTest)
            Output = p.map(funcSVR,param)
          
        #---------------------------------------------------------
        Output = np.asarray(Output)
        minArg = np.argmin(Output[:,0])
        #print(Output[minArg,1],Output[minArg,2],Output[minArg,3])
      
        gammaList.append(Output[minArg,1])
        cList.append(Output[minArg,2])
        epsilonList.append(Output[minArg,3])
        errorList.append(Output[minArg,0])
   
        filename2 = os.path.join('Experimento_3/data_save','Error','error'+str(gg)+'_'+str(ff)+'.npy')
        np.save(filename2,errorList)
    return cList, gammaList, epsilonList, errorList

###############################################################################
def Sample_Point_Cloud(filename,dv):
    pcdAux = opd.read_point_cloud(filename)
    auxPointsS = Feature_Scaling(np.asarray(pcdAux.points),-1,1)
    for pp,ii in enumerate(auxPointsS):
        pcdAux.points[pp] = ii
    pcd = opd.voxel_down_sample(pcdAux, voxel_size = dv)
    print('Point Cloud size =  ',np.asarray(pcd.dimension))
    
    pcdw = opd.PointCloud()
    pcdw.points = opd.Vector3dVector(pcd.points)
    opd.write_point_cloud('PointCloudSample.ply',pcdw)
    
##############################################################################
def get_data_myc(filename):
  
    #filename = os.path.join('Points_Clouds',name)
    pcd  = opd.read_point_cloud(filename)
    
    #auxPointsS = Feature_Scaling(np.asarray(pcd.points),-1,1)
    #for pp,ii in enumerate(auxPointsS):
    #    pcd.points[pp] = ii
      
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
###############################################################################

if __name__ == '__main__':
    #path = '/home/felipe/Desktop/Data_Heart/MM-WHS_2/Validacion_Surface/Point_cloud'
    nPart = 10
    filename = os.path.join('Point_Cloud','MycInferior.ply')
    Sample_Point_Cloud(filename,.09)
    pcdAux = opd.read_point_cloud('PointCloudSample.ply')
    
    #pcdAux  = opd.read_point_cloud(filename)
    Train_Index,Test_Index = crossval_val2( np.asarray(pcdAux.points),nPart)
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
        
        pointsAuxTrain = np.asarray(pcdAux.points)[Train_Index[ii]]
        pointsAuxTest = np.asarray(pcdAux.points)[Test_Index[ii]]
        
        pcd = opd.PointCloud()
        pcd.points = opd.Vector3dVector(pointsAuxTrain)
        opd.write_point_cloud('Temp.ply',pcd)
        
        Xtrain,tTrain = get_data_myc('Temp.ply')
        #-----------------------------------------------------

        cList, gammaList, epsilonList, errorList = FindParamSVR(
                                                    Xtrain,tTrain,Kernel,pointsAuxTrain,ii)
#        #--------------------------------------------
#        #Best Parametros
        indx = np.argmin(errorList)
#      
        cbest = cList[indx]
        gammabest = gammaList[indx]
        epsilonbest = epsilonList[indx]

        Xtest = pointsAuxTest 
        svr,t_est2 = svrTrain(Xtrain,tTrain,Xtest,Kernel,gammabest,cbest,epsilonbest)
        filename2 = os.path.join('Experimento_3/Model','svrMycSUP'+str(ii+1)+'.joblib')
        print('SAVE ',filename2,'\n')
        dump(svr, filename2)
        tTest2 = np.zeros((Xtest.shape[0],))

        ## Error RMSD
        Error2 =  np.sqrt((np.sum((t_est2-tTest2)**2)/tTest2.shape[0]))
      
        errorList2.append(Error2)
#      
#        #-------------------------------------------
#        #save
        cDic[str(ii)] = cList
        gammaDic[str(ii)] = gammaList
        epsilonDic[str(ii)] = epsilonList
        errorDic[str(ii)] = errorList
#      
        np.save('Experimento_3/data_save/cDic.npy',cDic)
        np.save('Experimento_3/data_save/gammaDic.npy',gammaDic)
        np.save('Experimento_3/data_save/epsilonDic.npy',epsilonDic)
        np.save('Experimento_3/data_save/errorDic.npy',errorDic)
        np.save('Experimento_3/data_save/error2Dic.npy',errorList2)
