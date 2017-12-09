'''
Created on Nov 13, 2017

@author: dduque
'''
from sklearn.cluster import KMeans

from builtins import int
import numpy as np
from ScenarioReduction.ScenarioObjects import Scenario
from multiprocessing import Pool

scenarios = None
J =None


def scenarioDistance(sc1, sc2):
        dif = sc1.scenarioVector - sc2.scenarioVector
        dist = np.linalg.norm(dif,1)
        return dist
    
def mp_solveFor_l(index):
    val = 0
    tempJ = dict(J)
    tempJ[index] = scenarios[index]
    for k in tempJ:
        sck = scenarios[k]
        pk = sck.p
        minIndexj = -1
        minVal = float("inf")
        gen_l = (x for x in range(0,len(scenarios)) if x not in tempJ)
        for j in gen_l:
            dist_kj = sck.computeDistance(scenarioDistance,scenarios[j])
            if dist_kj<minVal:
                minVal = dist_kj
                minIndexj = j
        val = val + pk*minVal
    return val, index

def mp_lstart(p, distFunc):
    global J
    global scenarios
    
    testValues = [s.id for s in scenarios if s.id not in J]
    myfunc = mp_solveFor_l
    ans = p.map( myfunc , testValues)
    minScenario = min(ans , key=lambda x: x[0])
    return minScenario[1]

def mp_computeReducedSet(scenarioList, n, pool_size = 10):
    reducedSet = {}
    global J
    J = {}
    global scenarios
    scenarios = scenarioList
    p = Pool(pool_size)
    steps = len(scenarios)-n
    for i in range(0,steps):
        l_start = mp_lstart(p , scenarioDistance)
        J[l_start] = scenarios[l_start]
    
    for sc in scenarios:
        if sc.id not in J:
            reducedSet[sc.id] = Scenario(sc.id,sc.scenarioVector, sc.p )
    
    for r in J:
        scr = J[r]
        targetSce = scr.findClosest(scenarioDistance, reducedSet)
        targetSce.p += scr.p
    return reducedSet
    
from sklearn import cluster as cl

    
def kmeans_reduce_set(scenarioList, n):
    #Build data set
    X = np.array([s.scenarioVector for s in scenarioList])
    resultKmean = cl.AgglomerativeClustering(n,'l1' ,linkage = 'average').fit(X)
    clusters = {k:{} for k in range(0,n)}
    for (i,s) in enumerate(scenarioList):
        clusters[resultKmean.labels_[i]][s.id] = s
    
    reducedSet = {}
    for c in clusters.values():
        sc = get_center_scenario(c)
        new_p = sum(s.p for s in c.values())
        reducedSet[sc.id] = Scenario(sc.id, sc.scenarioVector, new_p)
    return reducedSet
 
def get_center_scenario(c):
    norm_p = sum(s.p for s in c.values())
    cluster_vectors = np.array([s.p*s.scenarioVector/norm_p for s in c.values()])
    centroid_vals = np.mean(cluster_vectors, 0)
    fic_sce = Scenario(-1,centroid_vals,0)
    sce = fic_sce.findClosest(scenarioDistance, c)
    return sce   
    
    






class ScenarioReductionAlgo(object):
    '''
    classdocs
    '''


    def __init__(self, scenarioList = None , dFunction = None):
        '''
        Constructor
        '''
        if scenarioList != None:
            self.scenarios = scenarioList
        else:
            self.scenarios = []
            
        self.reducedSet = {}
        self.distFunction = dFunction 
    
    
    def computeReducedSet(self, n):
        '''
        Run a scenario reduction heuristic
        Args:
            n (int): target parameters
        '''
        J = {}
        steps = len(self.scenarios)-n+1
        for i in range(0,steps):
            print(i)
            l_start = self.getLstar(J, greedySearch = True)
            lstar2 = mp_lstart(J, self.scenarios, self.distFunction)
            J[l_start] = self.scenarios[l_start]
        
        for sc in self.scenarios:
            if sc.id not in J:
                self.reducedSet[sc.id] = Scenario(sc.id,sc.scenarioVector, sc.p )
        
        for r in J:
            scr = J[r]
            targetSce = scr.findClosest(self.distFunction, self.reducedSet)
            targetSce.p += scr.p
        
    
        return self.reducedSet
    def getLstar(self, J, greedySearch=False):
        '''
        Args:
            J (dict of Scenarios): Scenarios to remove
        '''
        minVal = float("inf")
        minIndex = -1
        for (l,sc) in enumerate(self.scenarios):
            if l not in J:
                val_l, index = self.getTestValue(l, J)
                if val_l<minVal:
                    minVal = val_l
                    minIndex = l
        return minIndex  
         
    def getLstar_mp(self, J, greedySearch=False):
        '''
        Args:
            J (dict of Scenarios): Scenarios to remove
        '''
        p = Pool(len(self.scenarios))
        gen_map = (l for (l,sc) in enumerate(self.scenarios) if l not in J)
        
        ans = p.map(getTestValJ(J,self), gen_map)
        minIndex =2
        return minIndex  
    
    def getTestValue(self, index, J):
        val = 0
        tempJ = dict(J)
        tempJ[index] = self.scenarios[index]
        for k in tempJ:
            sck = self.scenarios[k]
            pk = sck.p
            minIndexj = -1
            minVal = float("inf")
            gen_l = (x for x in range(0,len(self.scenarios)) if x not in tempJ)
            for j in gen_l:
                dist_kj = sck.computeDistance(self.distFunction,self.scenarios[j])
                if dist_kj<minVal:
                    minVal = dist_kj
                    minIndexj = j
            val = val + pk*minVal
        return val, index
                
            
                
            
            
        
        #Clean up
        del(tempJ)
       
    
