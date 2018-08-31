# -*- coding: utf-8 -*-
"""
Created on Thu Aug 09 10:47:36 2018

@author: 913
"""

import triangle
import triangle.plot as triplt
import numpy as np

import itertools
import matplotlib.pylab as plt
import pandas as pd



def saveTriaPoly(dictPoly,name):
    assert ("vertices" in dictPoly.keys())
    # 保存节点
    try:
        if "vertex_markers" in dictPoly.keys():
            lines = list()
            lines.append("%d 2 0 1\n"%(dictPoly['vertices'].shape[0]))
            for ii in xrange(dictPoly['vertices'].shape[0]):
                cell = dictPoly['vertices'][ii]
                lines.append("%d %.6f %.6f %d\n"%(ii+1, cell[0], cell[1], dictPoly['vertex_markers'][ii]))
        else:
            lines = list()
            lines.append("%d 2 0 0\n"%(dictPoly['vertices'].shape[0]))
            for ii in xrange(dictPoly['vertices'].shape[0]):
                cell = dictPoly['vertices'][ii]
                lines.append("%d %.6f %.6f %d\n"%(ii+1, cell[0], cell[1], 0))            
    except Exception as e:
        print "!"*8,'vertices<<--',e
        raise
    try: 
        # 保存线段
        if "segment_markers" in dictPoly.keys():           
            lines.append("%d 1\n"%dictPoly['segments'].shape[0])
            for ii in xrange(dictPoly['segments'].shape[0]):
                cell = dictPoly['segments'][ii]
                lines.append("%d %d %d %d\n"%(ii+1, cell[0]+1, cell[1]+1, dictPoly['segment_markers'][ii]))
            
        else:
            lines.append("%d 0\n"%dictPoly['segments'].shape[0])
            for ii in xrange(dictPoly['segments'].shape[0]):
                cell = dictPoly['segments'][ii]
                lines.append("%d %d %d\n"%(ii+1, cell[0]+1, cell[1]+1))
            
    except Exception as e:
        print "!"*8,'segment<<--',e
        lines.append('0')
        
    try:
        lines.append("%d\n"%dictPoly['holes'].shape[0])
        # 保存空洞位置
        for ii in xrange(dictPoly['holes'].shape[0]):
            cell = dictPoly['holes'][ii]
            lines.append("%d %.6f %.6f\n"%(ii+1, cell[0], cell[1]))
    except Exception as e:
        print "!"*8, 'holes<<--',e
        lines.append('0')
    try:
        # 文件操作
        with open(name+".poly",'w') as f:
            f.writelines(lines)
            pass
    except Exception as e:
        print "!"*8, 'vertices<<--',e