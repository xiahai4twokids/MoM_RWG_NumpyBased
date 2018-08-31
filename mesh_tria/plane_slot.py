#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 10:51:16 2018

@author: junhu
"""

import itertools
import triangle
import triangle.plot as triplt
import numpy as np
import matplotlib.pylab as plt
from SaveTria import *



# 生成节点
plan_width = 1.5
xgrid = [-plan_width,plan_width]
ygrid = [-plan_width,plan_width]
grid4corner = np.array(list(itertools.product(xgrid,ygrid)))

slot_length = 1
slot_width = 0.03
xgrid_slot = [-slot_length,slot_length]
ygrid_slot = [-slot_width,slot_width]
grid4slot1 = np.array(list(itertools.product(xgrid_slot,ygrid_slot)))
grid4slot2 = np.dot(grid4slot1,np.array([[0,1],[-1,0]]))

xgrid_inner_conner = ygrid_slot
ygrid_inner_conner = ygrid_slot
grid4slot_inner_conner = np.array(list(itertools.product(xgrid_inner_conner,ygrid_inner_conner)))

vertices = np.vstack([grid4corner,grid4slot1,grid4slot2,grid4slot_inner_conner])
# 生成线段
segment1 = np.array([[0,1],[1,3],[3,2],[2,0]])
segment2 = np.array([[0,1],[1,9],[9,7],[7,6],
                     [6,11],[11,3],[3,2],[2,10],
                     [10,4],[4,5],[5,8],[8,0]])+4
segment = np.vstack([segment1,segment2])

# 生成空洞
hole = np.array([0,0]).reshape([1,-1])
# 完成PLSG的结构体
myPloy = dict({"vertices":vertices,'holes':hole,'segments':segment})


# 画出PLSG
plt.figure()
ax1 = plt.subplot(111, aspect='equal')
triplt.plot(ax1,**myPloy)
plt.show()

# 进行网格剖分
plt.figure()
myMesh = triangle.triangulate(myPloy,'pq30a0.01')
ax1 = plt.subplot(111, aspect='equal')
triplt.plot(ax1,**myMesh)
plt.show()

# 将剖分网格进行保存
saveTriaPoly(myMesh,"plane_slot")

# 从文件中导入网格
C_copy = triangle.load('.','plane_slot')
C_copy = triangle.triangulate(C_copy,'p') # 由于PLSG文件中不包含三角形，故需要进行组合成三角形
ax1 = plt.subplot(111, aspect='equal')
triplt.plot(ax1,**C_copy)
plt.show()
