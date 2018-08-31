# -*- coding: utf-8 -*-
# common packages
import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
from scipy.sparse.linalg import LinearOperator
import datetime
import time
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import pandas as pds

# defined packages
from _domain_utils import autoDecomDomain, optDecomDomain,optDecomDomain_check
#import pyximport; pyximport.install()
from Components import FillingProcess,FillingProcess_aca
from Components import ImpMatrix2,ImpMatrix3,ImpMatrix1
from Components import Solver,getFarFiled
from Components import RWGFunc

# In[] Some common parameters
from Parameters import Filename,DomainSeg,SolverPar,RCSPar_theta, RCSPar_phi, WorkingFreqPar

# In[]
def plotGeo(grids,trias):
    try:
        assert trias.shape[1] >= 3        
        x = grids[:,0]
        y = grids[:,1]
        z = grids[:,2]
        
        if trias.shape[-1] == 4:
            domainID = np.unique(trias[:,3])
            domainID_attached = pds.Series(domainID)
            center_domain = np.zeros([len(domainID),3])
            num_domain = np.zeros([len(domainID),1])
            for tria in trias:
                index_tria = domainID_attached[domainID_attached==tria[3]].index.values[0]
                center_domain[index_tria] = center_domain[index_tria] + \
                    np.array([np.mean(x[tria[:-1]]),np.mean(y[tria[:-1]]),np.mean(z[tria[:-1]])])
                num_domain[index_tria] = num_domain[index_tria]+1
            center_domain  = center_domain/num_domain    
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
        colors_key = colors.keys()
        if trias.shape[-1] == 4:
            try:
                for tria in trias:
                    indx = np.hstack((tria[:-1],tria[0]))
                    ax.plot(x[indx], y[indx], z[indx],color=colors_key[0])
            except Exception as e:
                print e
                raise
            for ii in xrange(center_domain.shape[0]):
                ax.text(center_domain[ii,0],center_domain[ii,1],center_domain[ii,2], 'D-%d'%ii,color='black')
        else:
            try:
                for tria in trias:
                    indx = np.hstack((tria,tria[0]))
                    ax.plot(x[indx], y[indx], z[indx],color=colors_key[0])
            except Exception as e:
                print e
                raise
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.show()
        
        pass
    except AssertionError as ae:
        print ae
        raise
    except Exception as e:
        print e
        raise
        
def plotDomain(grids,trias):
    try:
        assert trias.shape[1] >= 3        
        x = grids[:,0]
        y = grids[:,1]
        z = grids[:,2]
        
        if trias.shape[-1] == 4:
            domainID = np.unique(trias[:,3])
            domainID_attached = pds.Series(domainID)
            center_domain = np.zeros([len(domainID),3])
            num_domain = np.zeros([len(domainID),1])
            for tria in trias:
                index_tria = domainID_attached[domainID_attached==tria[3]].index.values[0]
                center_domain[index_tria] = center_domain[index_tria] + \
                    np.array([np.mean(x[tria[:-1]]),np.mean(y[tria[:-1]]),np.mean(z[tria[:-1]])])
                num_domain[index_tria] = num_domain[index_tria]+1
            center_domain  = center_domain/num_domain
        
                    
        
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
        colors_key = colors.keys()
        if trias.shape[-1] == 4:
            try:
                for tria in trias:
                    ax.plot_trisurf(x[tria[:-1]], y[tria[:-1]], z[tria[:-1]], \
                                    color = colors_key[tria[-1]%len(colors_key)], linewidth=0.3, antialiased=True)
            except Exception as e:
                print e
                raise
        else:
            try:
                for tria in trias:
                    ax.plot_trisurf(x[tria], y[tria], z[tria], color = colors_key[0], linewidth=0.3, antialiased=True)
            except Exception as e:
                print e
                raise                
        
        ax.set_xlabel('x')
        ax.set_ylabel('z')
        ax.set_zlabel('y')
        plt.show()
        
        pass
    except AssertionError as ae:
        print ae
        raise
    except Exception as e:
        print e
        raise
                
def preCal(wavenumber, grids, trias, edges, segment):              
    try:     
        trias__,domainGrid, domainIDs  = autoDecomDomain(grids, trias, segment)
        
        if len(optDecomDomain_check(8, grids, trias__,domainGrid, domainIDs)) != 0:
            optDecomDomain(8, grids, trias__,domainGrid, domainIDs)
        
        print "===="*30
        trias = trias__
        domains = np.unique(np.array(trias)[:,3])
        
        
        # grouping triangles
        triasinDomain = [[id for id,tria in enumerate(trias) if tria[3]==domainNum] for domainNum in domains]
        gridinDomain = [] # recording all nodes of domains, for determining neighborhood
        for domainNum in domains:
            temp = set()
            for tria in trias:
                if tria[3] == domainNum:
                    temp.update(tria[:-1])
            gridinDomain.append(temp)
    except Exception as e:
        print e
        raise   
    try:  
        k = wavenumber # wavenumber
    except Exception as e:
        print e
        raise         
    try:
        # qualifying mesh
        edge_leng = [np.linalg.norm(np.array(grids[edge[0]])-np.array(grids[edge[1]])) for edge in edges] 
        max_len = np.max(np.array(edge_leng))
        if k*max_len < scipy.pi*2*0.2: print 'good mesh'
        else: print 'poor mesh' 
    except Exception as e:
        print 'Skip', e,'-->> Edge'
         
    try:
        ################################################################## 
        # generating HRWG and RWG
        rwgs = RWGFunc().match(grids, trias)
        return [k,grids,trias,rwgs,domains,gridinDomain,triasinDomain]
    except Exception as e:
        print e
        raise
        
# In[]
def loadMem(filenamePar):
    try:
        # load structure from file
        import triangle as triaPy
        temp = triaPy.load('./',filenamePar.filename)
        temp = triaPy.triangulate(temp,'p') # generating triangles from PLSG file
        grids = temp['vertices']
        
        trias = temp['triangles']
        edges = temp['segments']
        # change structure
        grids = np.hstack( (grids[:,0:1], np.zeros([grids.shape[0],1]), grids[:,1:2]) )

        return [grids,trias,edges]
    except Exception as e:
        print e
        raise
# In[]
def simulator(filename=Filename(),solverPar=SolverPar()):
    try:
        sim_start = datetime.datetime.now()
        ID = "%s_%dH_%dM_%dS_%s_%s"%(\
                                     sim_start.date(),\
                                     sim_start.hour,\
                                     sim_start.minute,\
                                     sim_start.second,\
                                     filename.filename,\
                                     solverPar.matrix_solver_type\
                                     )
        details = dict()
        
        print "load mesh"
        grids, trias, edges = loadMem(filename)
                
        wavenumber=WorkingFreqPar().wavenumber
        segment = DomainSeg().segment
        
        k,grids,trias,rwgs,domains,gridinDomain,triasinDomain = \
            preCal(wavenumber, grids, trias, edges,segment)   

        plotGeo(grids,trias)
        grids_temp = np.hstack( (grids[:,0:1], grids[:,2:3], np.zeros([grids.shape[0],1])) )
        plotDomain(grids_temp,trias)
        
        details['Triangles'] = trias.shape[0]
        details['RWG'] = len(rwgs[1])
        details['freq'] = 3e8*k/np.pi/2
        details['wavenumber'] = k
        details['domains'] = len(domains)
        details['triaInDomain'] = triasinDomain
        
        print "Geo-Info"
        print "Triangles  = %d "%(details['Triangles'])        
        print "RWG Funcs = %d "%(details['RWG'])
        print 'frequency = %.2e Hz'%(details['freq'])
        print "wavenumber = %.2e"%(details['wavenumber'])
        print "domains = %d"%details['domains']
        
    except Exception as e:
        print e
        raise
    try:
        matrix_solver_type = solverPar.matrix_solver_type
        details['Solver_Info'] = matrix_solver_type
        print "Solver_Info = %s"%details['Solver_Info']
                
        print "filling matrix"
        filling_start = datetime.datetime.now()
        filling_cpu_start = time.clock()
        print 'filling start @ ', filling_start
        try:
            if matrix_solver_type == 'dir':
                matrix_all_in, filling_hander = FillingProcess().fillingProcess(k,grids,trias,rwgs,domains,gridinDomain,triasinDomain)
                matrix = LinearOperator((len(matrix_all_in[1]),len(matrix_all_in[1])),
                                        ImpMatrix2(matrix_all_in).matVec,ImpMatrix2(matrix_all_in).rmatVec)
            else:
                matrix_all_in, filling_hander = FillingProcess_aca().fillingProcess_aca(k,grids,trias,rwgs,domains,gridinDomain,triasinDomain)
                matrix = LinearOperator((len(matrix_all_in[1]),len(matrix_all_in[1])),
                                        ImpMatrix3(matrix_all_in).matVec,ImpMatrix3(matrix_all_in).rmatVec)
        except Exception as e:
            print e
            print "=*="*30
            raise
        filling_end = datetime.datetime.now()
        filling_cpu_end = time.clock()
        details['fillingtime'] = (filling_end-filling_start).seconds
        details['fillingcputime'] = filling_cpu_end-filling_cpu_start
        print 'filling end @ ', filling_end
        print 'filling time = %.2e s = %.2e m' %(details['fillingtime'], details['fillingtime']/60.)
        print "cpu time = %.2e s = %.2e m"%(details['fillingcputime'], details['fillingcputime']/60.)
        print "--"*90
    except Exception as e:
        print e
        raise
    except AssertionError as ae:
        print ae
        raise   
        
    try:
        print "solving equation"
        solving_start = datetime.datetime.now()
        solving_cpu_start = time.clock()
        print 'solving start @ ', solving_start
        result1 = Solver().cgsSolve(matrix, matrix_all_in[1]) 
        solving_end = datetime.datetime.now()
        solving_cpu_end = time.clock()
        details['solveingtime'] = (solving_end-solving_start).seconds
        details['solveingcputime'] = solving_cpu_end-solving_cpu_start
        details['rhd'] = matrix_all_in[1]

        print 'solving end @ ', solving_end
        print 'solving time =  %.2e s = %.2e m'%(details['solveingtime'], details['solveingtime']/60.)
        print "cpu time = %.2e s = %.2e m"%(details['solveingcputime'], details['solveingcputime']/60.)

        assert  result1[1] == 0
        I_current = result1[0].reshape([-1,1])
        details['current'] = I_current
        print "--"*90
    except Exception as e:
        print e
        raise
    except AssertionError as ae:
        print ae
        raise
   
    try:     
        print "calculating field on H-plane"        
        tempRCSPar = RCSPar_theta()
        theta_0 = tempRCSPar.theta_0
        phi_0 = tempRCSPar.phi_0
        r = tempRCSPar.r
        r_obs = [[[r*np.sin(theta)*np.cos(phi),r*np.sin(theta)*np.sin(phi),r*np.cos(theta)] \
                   for theta in theta_0] \
                        for phi in phi_0]
        theta_vector = [[[np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),-np.sin(theta)] \
                          for theta in theta_0] \
                            for phi in phi_0]
        phi_vector = [[[-np.sin(phi),np.cos(phi),0.] \
                        for theta in theta_0] \
                            for phi in phi_0]
        r_obs = np.array(r_obs)
        theta_vector = np.array(theta_vector)
        phi_vector = np.array(phi_vector)
        field_obs = getFarFiled(r_obs, I_current, filling_hander, trias, rwgs)
        print "--"*90
    except AssertionError as ae:
        print ae
        raise
        
    try:        
        #1
        field_theta = np.multiply(field_obs,theta_vector)
        field_theta = np.sum(field_theta, axis=2)
        field_theta = np.multiply(field_theta,np.conj(field_theta))
        aug = np.abs(field_theta)*r**2*4*np.pi
        aug = np.log10(aug)*10
        print "theta_comp"
        import matplotlib.pylab as plt
        for ii in xrange(len(phi_0)):
            plt.plot(theta_0*180/np.pi, aug[ii,:],label='$\\phi$=%.1f'%(phi_0[ii]*180/np.pi))
        plt.legend()
        plt.grid()
        plt.ylim(ymin=-30)
        plt.xlabel('$\\theta$ Degree')
        plt.ylabel('10log10(.) dB')
        plt.show()  
        
        details["theta"] = theta_0
        details['f_theta'] = aug

    except AssertionError as ae:
        print ae
        raise     
        
    try:     
        print "calculating field on E-plane"        
        tempRCSPar = RCSPar_phi()
        theta_0 = tempRCSPar.theta_0
        phi_0 = tempRCSPar.phi_0
        r = tempRCSPar.r
        r_obs = [[[r*np.sin(theta)*np.cos(phi),r*np.sin(theta)*np.sin(phi),r*np.cos(theta)] \
                   for theta in theta_0] \
                        for phi in phi_0]
        theta_vector = [[[np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),-np.sin(theta)] \
                          for theta in theta_0] \
                            for phi in phi_0]
        phi_vector = [[[-np.sin(phi),np.cos(phi),0.] \
                        for theta in theta_0] \
                            for phi in phi_0]
        r_obs = np.array(r_obs)
        theta_vector = np.array(theta_vector)
        phi_vector = np.array(phi_vector)
        field_obs = getFarFiled(r_obs, I_current, filling_hander, trias, rwgs)
        print "--"*90
    except AssertionError as ae:
        print ae
        raise
        
    try:  
        field_theta = np.multiply(field_obs,theta_vector)
        field_theta = np.sum(field_theta, axis=2)
        field_theta = np.multiply(field_theta,np.conj(field_theta))
        aug = np.abs(field_theta)*r**2*4*np.pi
        aug = np.log10(aug)*10        
        #3
        print "theta_comp"
        import matplotlib.pylab as plt
        for ii in xrange(len(theta_0)):
            plt.plot(phi_0*180/np.pi, aug[:,ii],label='$\\theta$=%.1f'%(theta_0[ii]*180/np.pi))
        plt.legend()
        plt.grid()
        plt.xlabel('$\\phi$ Degree')
        plt.ylabel('10log10(.) dB')
        plt.show()    
        
        details["phi"] = phi_0
        details['f_phi'] = aug
        
    except AssertionError as ae:
        print ae
        raise
    return (ID,details)
        
if __name__ == '__main__':
    import os
    import pickle
    try:
        os.mkdir('result')
    except Exception as e:
        print e
        
    ID_sim,details_sim = simulator(filename=Filename('butterfly'),solverPar=SolverPar('aca'))
    with open('result/%s.txt'%ID_sim,'w') as f:
        pickle.dump(details_sim,f)

    

        
