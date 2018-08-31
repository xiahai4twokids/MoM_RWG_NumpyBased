# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 16:37:34 2018

@author: 913
"""

import itertools
import triangle
import triangle.plot as triplt
import numpy as np
import matplotlib.pylab as plt
from SaveTria import saveTriaPoly


def obtainFullTria():
    try:
        a = 1
        b = 1
        
        # 预定义的格子店
        x = np.array([0,-a,a,a,-a])
        y = np.array([0,b,b,-b,-b])
        
        grid4 = np.array([(xx,yy) for xx,yy in zip(x,y)])
        edges = np.array([[0,1],[1,2],[2,0],\
                          [0,3],[3,4],[4,0]])
        myPloy_with_grid = dict({'vertices':grid4,'segments':edges})
        
        # 画出PLSG
        plt.figure()
        ax1 = plt.subplot(221, aspect='equal')
        triplt.plot(ax1,**myPloy_with_grid)
        
        # 进行网格剖分
        myMesh_tria = triangle.triangulate(myPloy_with_grid,'pq30a0.01')
        ax2 = plt.subplot(222, aspect='equal')
        triplt.plot(ax2,**myMesh_tria)
        plt.show()
                
        return myMesh_tria
    except Exception as e:
        print e
        


myMesh_tria = obtainFullTria()

saveTriaPoly(myMesh_tria,"../butterfly")
