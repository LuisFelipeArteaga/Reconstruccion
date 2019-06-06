#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 14:50:16 2019

@author: felipe
"""

import numpy as np
import pcl
import os
#import cv2 
#import random
#import pptk
from pyntcloud import PyntCloud
#import pandas as pd
#import matplotlib.pyplot as plt
#import skimage.morphology as skm
#import scipy.ndimage as spn
from Surface_distance.metricsSurface import*
from hausdorff import hausdorff
#import pcl.visualization

def DownsamplingVoxelGrid():
    return 0

def FeatureScaling(X,a,b):
    return a + ((X-np.min(X))*(b-a))/(np.max(X)-np.min(X))

def IterativeCloserPoint(points_in,points_out):
    
    cloud_in = pcl.PointCloud()
    cloud_in.from_array(points_in)
    cloud_out = pcl.PointCloud()
    cloud_out.from_array(points_out)
    
    icp = cloud_in.make_IterativeClosestPoint()
    converged, transf, estimate, fitness = icp.icp(cloud_in, cloud_out)
    
    return converged, transf, estimate, fitness



if __name__=='__main__':
    
#    PointCloud1 = pcl.load('Point_Cloud/Myc.ply')
#    points1 = PointCloud1.to_array()
#    points1 = FeatureScaling(points1,-1,1)
#    PointCloud1.from_array(points1)
#    
#    sor = PointCloud1.make_voxel_grid_filter()
#    sor.set_leaf_size(0.02,0.02,0.02)
#    cloud_filtered = sor.filter()
#    
#    
#    pcl.save(cloud_filtered,'cloud_filtered.ply')
    
    
#    PointCloud2 = pcl.load('Surface_Myc/svrMyc1.ply')
#    points2 = PointCloud2.to_array()
    
    
    
    pointsGrid1 = PyntCloud.from_file('Point_Cloud/Myc.ply')
    Points1 = pointsGrid1.points
    X1 = FeatureScaling(Points1.values,-1,1)
    Points1['x'] = X1[:,0];Points1['y'] = X1[:,1];Points1['z'] = X1[:,2]
    pointsGrid1.points = Points1   
    del X1, Points1
    
    Distance_Average = []
    Distance_Hausdorff1 = []
    Distance_Hausdorff2 = []

    path1 = 'Surface_Experimento_2_Inferior'
    for ii in range(10):
        pointsGrid2 = PyntCloud.from_file(os.path.join(path1,'svrMyc'+str(ii+1)+'.ply'))
        Points2 = pointsGrid2.points
        X2 = FeatureScaling(Points2.values,-1,1)
        Points2['x'] = X2[:,0];Points2['y'] = X2[:,1];Points2['z'] = X2[:,2]
        pointsGrid2.points = Points2
        del X2, Points2 
        
        
        #pointsGrid2.plot(mesh=False, backend="threejs")
        voxelgrid_id1 = pointsGrid1.add_structure("voxelgrid", n_x=100, n_y=100, n_z=100)
        voxelgrid_id2 = pointsGrid2.add_structure("voxelgrid", n_x=100, n_y=100, n_z=100)
        voxelgrid1 = pointsGrid1.structures[voxelgrid_id1]
        voxelgrid2 = pointsGrid2.structures[voxelgrid_id2]
        voxelgrid1.plot(d=3, mode="density", cmap="cool")
        voxelgrid2.plot(d=3, mode="density", cmap="cool")
        
    #    
        mask_gt = np.zeros((100, 100, 100), np.uint8)
        mask_pred = np.zeros((100, 100, 100), np.uint8)
        mask_gt[voxelgrid1.voxel_x,voxelgrid1.voxel_y,voxelgrid1.voxel_z] = 1
        mask_pred[voxelgrid2.voxel_x,voxelgrid2.voxel_y,voxelgrid2.voxel_z] = 1
        
        #Average Surface Distance
        surface_distances = compute_surface_distances(mask_gt, mask_pred, spacing_mm=(2,2,2))
        Average_Surface = compute_average_surface_distance(surface_distances)
    
        # Hausdorff Distance 
        Hausdorff_distance1 = compute_robust_hausdorff(surface_distances,100.0)
        Hausdorff_distance2 = hausdorff(pointsGrid1.points.to_numpy().astype(np.double),
                                       pointsGrid2.points.to_numpy().astype(np.double),
                                       distance = 'euclidean' )
        
        Distance_Average.append(Average_Surface)
        Distance_Hausdorff1.append(Hausdorff_distance1)     
        Distance_Hausdorff2.append(Hausdorff_distance2)
        
        np.save('Distance_Average'+path1[22:]+'.npy',Distance_Average)
        np.save('Distance_Hausdorff1'+path1[22:]+'.npy',Distance_Hausdorff1)
        np.save('Distance_Hausdorff2'+path1[22:]+'.npy',Distance_Hausdorff2)
        
        print('\n\n')
        print('Average Surface Distance ' , Average_Surface)
        print('Hausdorff Distance 1 ' , Hausdorff_distance1)
        print('Hausdorff Distance 2 ' , Hausdorff_distance2)
        #####################################################################################
        
        
    #    struct1 = ndimage.generate_binary_structure(3,2)
    #    struct2 = ndimage.generate_binary_structure(3,2)
    #    mask_gt = ndimage.binary_dilation(mask_gt,structure=struct1).astype(mask_gt.dtype)
    #    mask_gt = ndimage.binary_erosion(mask_gt,structure=struct2).astype(mask_gt.dtype)
    #    
    #    ab = mask_gt[:,:,51]
    #    kernel =  np.ones((3,3),np.uint8)
    #    kernel1 =  np.ones((10,10),np.uint8)
    #    struct3 = np.zeros((5,5),np.bool)
    #    ab2 = cv2.dilate(ab,kernel)
    #    ab3 = cv2.morphologyEx(ab2, cv2.MORPH_CLOSE, kernel1)
    #    ab4 = ndimage.binary_closing(ab2,structure=struct3).astype(ab.dtype)
        
    #    mask_gt_b = skm.binary_dilation(mask_gt,skm.cube(10))
    #    mask_gt_b2 = np.zeros((100, 100, 100), np.uint8)
    #    for ii in range(100):
    #        mask_gt_b2[:,:,ii] = spn.binary_fill_holes(mask_gt_b[:,:,ii]).astype(int)
    #        plt.figure(1)
    #        plt.imshow(mask_gt_b[ii,:,:])
    #        plt.show()
    
        
        
        
        
        
        




