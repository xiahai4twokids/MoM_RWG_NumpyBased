# -*- coding: utf-8 -*-


import numpy as np
import itertools
import pandas as pds

def autoDecomDomain(grids,trias, segment):
    try:
        if np.min(np.array([segment[0],segment[1], segment[2]])) > 0 :
            mins_grids = [np.min(np.array(grids)[:,k]) for k in xrange(3)] 
            maxs_grids = [np.max(np.array(grids)[:,k]) for k in xrange(3)] 
            start = np.array(mins_grids)
            end = np.array(maxs_grids)
            for xx in xrange(3):
                if start[xx] == end[xx]:
                    start[xx]=start[xx]-segment[xx]*1.e-8
                    end[xx] = end[xx]+segment[xx]*1.e-8
            # 计算三角形的中心
            center_cal = lambda tria : [np.mean([grids[tria[0]][k], grids[tria[1]][k], grids[tria[2]][k]]) for k in xrange(3)]
            centers = map(center_cal, trias)
    #        # 计算三个方向的极值坐标
    #        mins = [np.min(np.array(centers)[:,k]) for k in xrange(3)]
    #        maxs = [np.max(np.array(centers)[:,k]) for k in xrange(3)]
            # 计算三个方向的分段办法
            NofSeg = [int(np.ceil((end[ii]-start[ii])/segment[ii])) for ii in xrange(3)]
            intervals_1d = lambda ii : [(start[ii]+k*segment[ii], start[ii]+(k+1)*segment[ii]) for k in xrange(NofSeg[ii])]
            intervals_3d = map(intervals_1d,xrange(3)) 
            # 将三个方向的分区进行笛卡尔积
            intervalsID = [ tuple(xx) for xx in itertools.product(xrange(NofSeg[0]),xrange(NofSeg[1]),xrange(NofSeg[2]))]
            intervals_ = {xx:ii for ii,xx in enumerate(intervalsID)} # 建立字典结构的备份，便于HASH快速查找编号 
            # 根据三角形中心位置判断三角形所属区域
            downInt_1d = lambda x: int(np.floor((x[0]-x[1])/x[2])) # 计算区间下限的函数
            downInt_3d = lambda grid: map(downInt_1d,zip(grid,start,segment)) # 计算立方体块三个下限值的函数
            centers_downInt = map(downInt_3d, centers) # 计算三角形中心所属的区域下限值
            centers_domainID = [ intervals_[tuple(kk)] for kk in centers_downInt ]  # 利用分区的Hash结构判断所属区域的编号
            
            temp = np.array(trias)
            trias_domain = np.array(zip(temp[:,0],temp[:,1],temp[:,2],centers_domainID) )# 将三角形的区域编号代替原来的三角形属性值
            return [trias_domain, intervals_3d, intervalsID]
        else:
            mins_grids = [np.min(np.array(grids)[:,k]) for k in xrange(3)] 
            maxs_grids = [np.max(np.array(grids)[:,k]) for k in xrange(3)] 
            start = np.array(mins_grids)
            end = np.array(maxs_grids)
#            for xx in xrange(3):
#                if start[xx] == end[xx]:
#                    start[xx]=start[xx]-segment[xx]*1.e-8
#                    end[xx] = end[xx]+segment[xx]*1.e-8
            # 计算三角形的中心
            center_cal = lambda tria : [np.mean([grids[tria[0]][k], grids[tria[1]][k], grids[tria[2]][k]]) for k in xrange(3)]
            centers = map(center_cal, trias)
#            # 计算三个方向的极值坐标
#            mins = [np.min(np.array(centers)[:,k]) for k in xrange(3)]
#            maxs = [np.max(np.array(centers)[:,k]) for k in xrange(3)]
            # 计算三个方向的分段办法
            NofSeg = [1. for _ in xrange(3)]
            intervals_1d = lambda ii : [(start[ii], end[ii]),]
            intervals_3d = map(intervals_1d,xrange(3)) 
            # 将三个方向的分区进行笛卡尔积
            intervalsID = [(0,0,0),]
#            intervalsID = [ tuple(xx) for xx in itertools.product(xrange(NofSeg[0]),xrange(NofSeg[1]),xrange(NofSeg[2]))]
#            intervals_ = {xx:ii for ii,xx in enumerate(intervalsID)} # 建立字典结构的备份，便于HASH快速查找编号 
            
            temp = np.array(trias)
            zeros__ = np.zeros_like(temp[:,0])
            trias_domain = np.array(zip(temp[:,0],temp[:,1],temp[:,2],zeros__) )# 将三角形的区域编号代替原来的三角形属性值
            return [trias_domain, intervals_3d, intervalsID]           
        pass
    except:
        raise

def optDecomDomain(min_tria, grids,trias_domain, intervals_3d, intervalsID):
    try:
        Domains = np.unique(trias_domain[:,3])  # 所有区域的全局编号
        Domains_attached = pds.Series(trias_domain[:,3]) 
        
        count_Trias = lambda ii:  [Domains[ii], Domains_attached[Domains_attached==Domains[ii]].count() ]
        numberTriainDomain = map(count_Trias, xrange(Domains.shape[0])) # 统计区域内三角形个数
        numberTriainDomain = pds.DataFrame(numberTriainDomain) 
        sat_Trias = lambda ii:  list(Domains_attached[Domains_attached==Domains[ii]].index.values)
        triaIDinDomain_local = map(sat_Trias, xrange(Domains.shape[0])) #统计区域内三角形的标号           

        select_cond = numberTriainDomain[1]<=min_tria
#        select_DomainID_global = numberTriainDomain[select_cond][0].values # 统计被遴选的区域编号
        select_DomainID_local = numberTriainDomain[select_cond].index.values # 统计被遴选的区域编号
        
        un_select_DomainID_global = numberTriainDomain[~(select_cond)][0].values # 统计没有被遴选的区域编号
        un_select_DomainID_local = numberTriainDomain[~(select_cond)].index.values # 统计没有被遴选的区域编号


        select_tria = lambda ii: triaIDinDomain_local[select_DomainID_local[ii]]
        select_Tria_ID = map(select_tria, xrange(len(select_DomainID_local))) # 统计所有被遴选的三角形编号
        id_tria_mod =  np.hstack(select_Tria_ID)
        
        #统计所有未被遴选的区域的中心
        intervalsID_copy =  np.array(intervalsID)

        intervalsID_unselect = intervalsID_copy[un_select_DomainID_local]    
        get_center_x = lambda id_iterval:  np.mean(intervals_3d[0][intervalsID_unselect[id_iterval][0]] )
        get_center_y = lambda id_iterval:  np.mean(intervals_3d[1][intervalsID_unselect[id_iterval][1]] )
        get_center_z = lambda id_iterval:  np.mean(intervals_3d[2][intervalsID_unselect[id_iterval][2]] )
        center_x = map(get_center_x, xrange(len(intervalsID_unselect)))
        center_y = map(get_center_y, xrange(len(intervalsID_unselect)))
        center_z = map(get_center_z, xrange(len(intervalsID_unselect)))
        centers_domains = np.array(zip(center_x,center_y,center_z))

        grid_tria = grids[trias_domain[id_tria_mod][:,0:3]]
        center_tria = np.mean(grid_tria,axis=1)
        
        distance_vec = center_tria.reshape([-1,1,3])-centers_domains.reshape([1,-1,3])
        distance = np.sum(distance_vec**2,axis=-1)
        targe_ID_domain_local =  np.argmin(distance,axis=1)

        targe_ID_domain_global = un_select_DomainID_global[targe_ID_domain_local]    
        trias_domain[id_tria_mod,3] = targe_ID_domain_global        
        Domains = np.unique(trias_domain[:,3])        
    except Exception as e:
        print e
        
        raise

def optDecomDomain_check(min_tria, grids,trias_domain, intervals_3d, intervalsID):
    try:
        Domains = np.unique(trias_domain[:,3])  # 所有区域的全局编号
        Domains_attached = pds.Series(trias_domain[:,3]) 
        
        count_Trias = lambda ii:  [Domains[ii], Domains_attached[Domains_attached==Domains[ii]].count() ]
        numberTriainDomain = map(count_Trias, xrange(Domains.shape[0])) # 统计区域内三角形个数
        numberTriainDomain = pds.DataFrame(numberTriainDomain) 
        
        select_cond = numberTriainDomain[1]<=min_tria
        select_DomainID_global = numberTriainDomain[select_cond][0].values # 统计被遴选的区域编号
        
        return select_DomainID_global
    except Exception as e:
        print e
        
        raise    
        
if __name__ == '__main__':
    import re
    # 读取网格文件
    f = open('B.1.node','r') # 读取网格节点文件
    head = f.readline()
    numbers = re.compile('\d+\.*\d*')
    NumofNodes = int(numbers.findall(head)[0])
    grids = []
    for str in f.readlines():
        m = numbers.findall(str)
   	if len(m)>2:
             grids.append((float(m[1]),float(m[2]), 0.))
    f.close() 
    g = open('B.1.ele','r') # 读取网格三角形文件
    head = g.readline()
    NumofTria = int(numbers.findall(head)[0])
    trias = []
    for str in g.readlines():
        m = numbers.findall(str)
        if len(m)>3:
            trias.append((int(m[1])-1,int(m[2])-1,int(m[3])-1, int(m[4])))
    g.close() 

    segment = [0.005,0.005,0.005]
    trias__,domainGrid = autoDecomDomain(grids, trias, segment)
    print np.array(trias__)
    print domainGrid
