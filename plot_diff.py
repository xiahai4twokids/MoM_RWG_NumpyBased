# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 14:31:02 2018

@author: 913
"""
# In[]
import os
try:
    os.mkdir('result')
except Exception as e:
    print e
try: 
    os.mkdir('tempData')
except Exception as e:
    print e

# In[]

import pickle
#from multiprocessing import Pool as ProgPool

from mom_solver import Solution
from mom_solver import Parameters

name = "butterfly"

#ID_sim_aca,details_sim_aca = Solution.simulator(filename=Parameters.Filename(name),solverPar=Parameters.SolverPar('aca'))
#with open('result/%s.txt'%ID_sim_aca,'w') as f:
#    print ID_sim_aca
#    pickle.dump(details_sim_aca,f)

ID_sim_dir,details_sim_dir = Solution.simulator(filename=Parameters.Filename(name),solverPar=Parameters.SolverPar('dir'))
with open('result/%s.txt'%ID_sim_dir,'w') as f:
    print ID_sim_dir
    pickle.dump(details_sim_dir,f)
    
#pool = ThreadPool(2)
#results = pool.map(testmulf,[1,2])
#pool.close()
#pool.join()
    
# In[]
#import matplotlib.pylab as plt
#plt.figure()
#plt.plot(details_sim_aca['theta'], details_sim_aca['f_theta'][0],'r-',label="aca")
#plt.plot(details_sim_dir['theta'], details_sim_dir['f_theta'][0],'k-',label='direct')
#plt.ylim(ymin=-20)
#plt.legend()
#plt.grid()
#plt.show()
#
#plt.figure()
#plt.plot(details_sim_aca['phi'], details_sim_aca['f_phi'][:,0],'r-',label="aca")
#plt.plot(details_sim_dir['phi'], details_sim_dir['f_phi'][:,0],'k-',label='direct')
#plt.ylim(ymin=-20)
#plt.legend()
#plt.grid()
#plt.show()

# In[]
#print details_sim_aca['fillingcputime']
#print details_sim_dir['fillingcputime']
#
#print details_sim_aca['solveingcputime']
#print details_sim_dir['solveingcputime']
