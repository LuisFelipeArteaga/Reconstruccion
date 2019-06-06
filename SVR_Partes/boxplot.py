#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 10:02:52 2019

@author: felipe
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



if __name__ == '__main__':
    path1 = 'Experimento_1_Inferior/data_save/'
    c= np.load((path1+'cDic.npy'))
    gamma = np.load((path1+'gammaDic.npy'))
    epsilon = np.load((path1+'epsilonDic.npy'))
    error = np.load((path1+'errorDic.npy'))
    error2 = np.load((path1+'error2Dic.npy'))
    
    d2f = 10
    auxBoxC = np.zeros((10,d2f))
    
    for jj in range(d2f):
        cAux = c.item().get(str(jj))
        auxBoxC[:,jj] = cAux
        print('Fold C= ',jj+1)
        print('Mean ',np.mean(cAux))
        print('Median ',np.median(cAux))
        print('std ',np.std(cAux),'\n')
    #print(cAux)
    #print(auxBoxC)
    plt.figure('C')
    df = pd.DataFrame(auxBoxC,
					columns =['c-svr1','c-svr2','c-svr3','c-svr4','c-svr5',
                                  'c-svr6','c-svr7','c-svr8','c-svr9','c-svr10'])
    boxplot = df.boxplot(column =['c-svr1','c-svr2','c-svr3','c-svr4','c-svr5',
                                  'c-svr6','c-svr7','c-svr8','c-svr9','c-svr10'])
    
    auxBoxgamma = np.zeros((10,d2f))
    for jj in range(d2f):
        cAux = gamma.item().get(str(jj))
        auxBoxgamma[:,jj] = cAux
        print('Fold Gamma= ',jj+1)
        print('Mean ',np.mean(cAux))
        print('Median ',np.median(cAux))
        print('std ',np.std(cAux),'\n')
    #print(cAux)
    #print(auxBoxgamma)
    plt.figure('Gamma')
    df = pd.DataFrame(auxBoxgamma,
					columns =['G-svr1','G-svr2','G-svr3','G-svr4','G-svr5',
                                  'G-svr6','G-svr7','G-svr8','G-svr9','G-svr10'])
    boxplot = df.boxplot(column =['G-svr1','G-svr2','G-svr3','G-svr4','G-svr5',
                                  'G-svr6','G-svr7','G-svr8','G-svr9','G-svr10'])
    
    auxBoxepsilon = np.zeros((10,d2f))
    for jj in range(d2f):
        cAux = epsilon.item().get(str(jj))
        auxBoxepsilon[:,jj] = cAux
        print('Fold Epsilon = ',jj+1)
        print('Mean ',np.mean(cAux))
        print('Median ',np.median(cAux))
        print('std ',np.std(cAux),'\n')
    #print(cAux)
    #print(auxBoxgamma)
    plt.figure('Epsilon')
    df = pd.DataFrame(auxBoxepsilon,
					columns =['E-svr1','E-svr2','E-svr3','E-svr4','E-svr5',
                                  'E-svr6','E-svr7','E-svr8','E-svr9','E-svr10'])
    boxplot = df.boxplot(column =['E-svr1','E-svr2','E-svr3','E-svr4','E-svr5',
                                  'E-svr6','E-svr7','E-svr8','E-svr9','E-svr10'])
    
    auxBoxerror = np.zeros((10,d2f))
    for jj in range(d2f):
        cAux = error.item().get(str(jj))
        auxBoxerror[:,jj] = cAux
        print('Fold Epsilon = ',jj+1)
        print('Mean ',np.mean(cAux))
        print('Median ',np.median(cAux))
        print('std ',np.std(cAux),'\n')
    #print(cAux)
    #print(auxBoxgamma)
    plt.figure('Error')
    df = pd.DataFrame(auxBoxerror,
					columns =['E-svr1','E-svr2','E-svr3','E-svr4','E-svr5',
                                  'E-svr6','E-svr7','E-svr8','E-svr9','E-svr10'])
    boxplot = df.boxplot(column =['E-svr1','E-svr2','E-svr3','E-svr4','E-svr5',
                                  'E-svr6','E-svr7','E-svr8','E-svr9','E-svr10'])