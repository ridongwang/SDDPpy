'''
Created on Oct 26, 2018

@author: Daniel Duque
'''
import numpy as np 

#CONSTANTS
demand_log_mean = 1.2
demand_log_std = 1.0

class Datahandler:
    '''
    This class manage the input data of the Capacity Allocation model 
    '''
    def __init__(self, n,m,scen,b,rho,seed = 0, dh=None):
        '''
        Define attributes used throughout the models
        '''
        #Sets
        self.Omega = None   #Set of scenarios
        self.I = None       #Set of generators
        self.J = None       #Set of Customers
        
        #Parameters
        self.prob = None    #Prob. of scenario w
        self.c  = None      #Unit operation cost of energy sent form generator i to customer j
        self.rho = None     #Unit subcontracting cost for demand site j
        self.b = None     #Maximum capacity instalation (scalar)
        
        '''Random parameters'''
        self.d = None       #Demand for customer j in scenario w
        
        if dh == None:
            self.buildData(n,m,scen,b,rho,seed)
        else:
            self.copy_dh(dh)
            self.generate_more_scenarios(scen, seed)
    
    def buildData(self, n,m,scen,b,rho, seed ):
        '''
        Builds an instance of the problem
        ==================
        input parameters:
        n: 
            Number of customers
        m: 
            Number of facilities
        scen:
            Number of scenarios
        b:
            butget
        rho: 
            Shortfall penalty
            
        #Original code in AMPL capacity.run
        #Author: Prof. David Morton
        '''
        self.Omega = list(range(0,scen))
        self.I = list(range(0,m))
        self.J = list(range(0,n))
        self.b = b
        self.rho = rho
        # pmf
        self.prob = np.ones(scen)*(1.0/scen)
        #Seed for locations
        np.random.seed(0)
        # Facility locations (x,y) coordinates in unit square
        xi = np.zeros(m)
        yi = np.zeros(m)
        for i in self.I:
            xi[i] = np.random.uniform()
            yi[i] = np.random.uniform()
        
        # Customer locations (x,y) coordinates in unit square
        xj = np.zeros(n)
        yj = np.zeros(n)
        for i in self.J:
            xj[i] = np.random.uniform()
            yj[i] = np.random.uniform()
        
        # unit cost of satisfying demand j from i is proportional to distance
        self.c = np.zeros((m,n))
        for i in self.I:
            for j in self.J:
                self.c[i,j]=np.round(np.sqrt((xi[i]-xj[j])**2+(yi[i]-yj[j])**2),3)
        
        #Set a seed
        np.random.seed(seed)
        # demands are independent lognormals
        self.d = np.zeros((n,scen))
        for w in self.Omega:
            for j in self.J:
                self.d[j,w] = np.round(np.exp(np.random.normal(demand_log_mean,demand_log_std)) , 3)
                
    def copy_dh(self, dh):
        self.__dict__ = dh.__dict__.copy()
        self.d = None 
    
    def generate_more_scenarios(self, scen, seed):
        #Set a seed
        np.random.seed(seed)
        self.Omega = list(range(0,scen))
        self.d = 0
        # demands are independent lognormals
        self.d = np.zeros((len(self.J),scen))
        for w in self.Omega:
            for j in self.J:
                self.d[j,w] = np.round(np.exp(np.random.normal(demand_log_mean,demand_log_std)) , 3)