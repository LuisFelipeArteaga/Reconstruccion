"""
Created on Fri Mar 15 14:50:16 2019

@author: felipe
"""

import numpy as np
import pcl
import pptk
from pyntcloud import PyntCloud
#import pcl.visualization

def DownsamplingVoxelGrid():
    return 0

def FeatureScaling(X,a,b):
    return a + ((X-np.min(X))*(b-a))/(np.max(X)-np.min(X))

if __name__=='__main__':
    
    PointCloud1 = pcl.load('Point_Cloud/Myc.ply')
    points1 = PointCloud1.to_array()
    
    
    PointCloud2 = pcl.load('Surface_Myc/svrMyc1.ply')
    points2 = PointCloud2.to_array()
    
    
    aa = PyntCloud.from_file('Point_Cloud/Myc.ply')
    
    #rgb  = pptk.rand(100,3)
    #v = pptk.viewer(points1)
    #v.set(point_size = 0.05)
    #v.wait()
