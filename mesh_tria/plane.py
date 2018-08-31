# -*- coding: utf-8 -*-
"""
Created on Thu Aug 09 10:46:44 2018

@author: 913
"""


import itertools
import triangle
import triangle.plot as triplt
import numpy as np
import matplotlib.pylab as plt
from SaveTria import saveTriaPoly
import pandas as pd

def obtainFullTria(x,y):
    try:
        y_temp,x_temp = np.meshgrid(y,x)
        grid4 = np.array(zip(x_temp.ravel(),y_temp.ravel()))
#        grid4 = np.array([(xx,yy) for xx in x for yy in y])
        myPloy_with_grid = dict({'vertices':grid4})
#        raise
        
        # 画出PLSG
        plt.figure()
        ax1 = plt.subplot(221, aspect='equal')
        triplt.plot(ax1,**myPloy_with_grid)
#        plt.show()
        
        # 进行网格剖分
        myMesh_tria = triangle.triangulate(myPloy_with_grid)
        ax2 = plt.subplot(222, aspect='equal')
        triplt.plot(ax2,**myMesh_tria)
#        plt.show()
        
        d = list(itertools.permutations([0,1,2],r=2))
        e = [myMesh_tria['triangles'][:,d[ii]] for ii in xrange(len(d))]
        edges = np.vstack(e)
        diff = edges[:,0]-edges[:,1]
        edges_2 = edges[diff>0,:]
        
        
        frame=pd.DataFrame(edges_2)
#        print frame
#        print  frame[frame.duplicated()]
        dup = list(frame.duplicated()==False)
        dup = [ii for ii,_ in enumerate(dup)]
        edges_3 = edges_2[dup,:]
#        print edges_3
#        print dup
        
        myPloy_with_seg = myPloy_with_grid
        myPloy_with_seg['segments'] = edges_3
        ax3 = plt.subplot(223, aspect='equal')
        triplt.plot(ax3,**myPloy_with_seg)
        for ii in xrange(myPloy_with_seg['vertices'].shape[0]):
            ax1.text(myPloy_with_seg['vertices'][ii,0],myPloy_with_seg['vertices'][ii,1],ii)
#        plt.show()
#        
        # 形成最终完整的mesh进行保存
        myMesh_tria_2 = triangle.triangulate(myPloy_with_seg,'p')
        ax4 = plt.subplot(224, aspect='equal')
        triplt.plot(ax4,**myMesh_tria_2)
        plt.show()
        
        return myMesh_tria_2
    except Exception as e:
        print e
        print grid4.shape
        raise

a = 2
b = 2

# 预定义的格子店
x = np.linspace(-a,a,20)
y = np.linspace(-b,b,20)

myMesh_tria = obtainFullTria(x,y)

saveTriaPoly(myMesh_tria,"../plane")
