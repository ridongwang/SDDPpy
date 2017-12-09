'''
Created on Nov 13, 2017

@author: dduque
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import matplotlib
font = {
        'weight' : 'bold',
        'size'   : 8}


matplotlib.rc('font', **font)

from bokeh.util.session_id import random
from ScenarioReduction.ScenarioObjects import Scenario, ScenarioSet
np.random.seed(0)
from ScenarioReduction.SRAlg import ScenarioReductionAlgo, mp_computeReducedSet, kmeans_reduce_set

def defineScenarioFunction(r): 
    '''
    Defines a distance
    '''  
    def scenarioDistance(sc1, sc2):
        dif = sc1.scenarioVector - sc2.scenarioVector
        dist = np.linalg.norm(dif, r)
        return dist
    
    return scenarioDistance

def sampledKdimensionalScenarios(k=1, N=100, covar = -1.5, mu = 1, var=4):
    '''
    Sample N k-dimensional vectors following a normal distribution.
    '''
    sigma = np.full((k,k), covar)
    for i in range(0,k):
        sigma[i,i] = var
    mu_vector = mu*np.ones(k)
    Z =  np.random.multivariate_normal(mu_vector,sigma,N)
    scenarioList = []
    for i in range(0,N):
        scenarioList.append(Scenario(i, Z[i], 1.0/N))
    return scenarioList , Z
    

if __name__ == '__main__':
    k = 2
    N = 500
    r = 1.0
    
    distFunction = defineScenarioFunction(r)
    scenarioList, randNums = sampledKdimensionalScenarios(k=k, N=N)
    meanVector = np.average(randNums, 0)
    covMatrix = np.cov(randNums.transpose())
    print(meanVector)
    print(covMatrix)
    f, axarr = plt.subplots(4,2,figsize=(5, 10), dpi=300)
    #gs1 = gridspec.GridSpec(4, 2)
    #gs1.update(wspace=0.9, hspace=10.5) # set the spacing between axes. 
    for (i,n)in enumerate([10,20,50,100]):
        print('Runun kmeans with n=%i' %(n))
        reducedset = kmeans_reduce_set(scenarioList, n)
        #sr_algorithm = ScenarioReductionAlgo(scenarioList, distFunction)
        #reducedset = mp_computeReducedSet(scenarioList, n)
        rSet = ScenarioSet(k, reducedset)
        print(rSet.checkProbConsisntency())
        print(rSet.computeExpectation())
        print(rSet.computeCovarianve())
        axij = axarr[i,1]
        x1 = [s.scenarioVector[0] for s in scenarioList]
        y1 = [s.scenarioVector[1] for s in scenarioList]
        s1 = [s.p*N for s in scenarioList]
        axij.scatter(x1,y1,s=s1, color='r')
        x2 = [s.scenarioVector[0] for s in reducedset.values()]
        y2 = [s.scenarioVector[1] for s in reducedset.values()]
        s2 = [s.p*N for s in reducedset.values()]
        axij.scatter(x2,y2,s=s2, color='b')
        axij.set_xlim([-5, 5])
        axij.set_ylim([-5, 5])
        axij.set_aspect('equal', 'box-forced')
        axij.set_title('k-means algorithm n=%i' %(n))
    
    for (i,n)in enumerate([10,20,50,100]):
        print('Runing Backward reduction with n=%i' %(n))
        #reducedset = kmeans_reduce_set(scenarioList, n)
        #sr_algorithm = ScenarioReductionAlgo(scenarioList, distFunction)
        reducedset = mp_computeReducedSet(scenarioList, n)
        rSet = ScenarioSet(k, reducedset)
        print(rSet.checkProbConsisntency())
        print(rSet.computeExpectation())
        print(rSet.computeCovarianve())
        axij = axarr[i,0]
        x1 = [s.scenarioVector[0] for s in scenarioList]
        y1 = [s.scenarioVector[1] for s in scenarioList]
        s1 = [s.p*N for s in scenarioList]
        axij.scatter(x1,y1,s=s1, color='r')
        x2 = [s.scenarioVector[0] for s in reducedset.values()]
        y2 = [s.scenarioVector[1] for s in reducedset.values()]
        s2 = [s.p*N for s in reducedset.values()]
        axij.scatter(x2,y2,s=s2, color='b')
        axij.set_xlim([-5, 5])
        axij.set_ylim([-5, 5])
        axij.set_aspect('equal', 'box-forced')
        axij.set_title('Backward algorithm n=%i' %(n))
    #f.tight_layout()
    plt.draw()
    plt.subplots_adjust(wspace=0.5, hspace=0.4)
    pp = PdfPages('./sr_figs.pdf')
    pp.savefig(f)
    pp.close()
    plt.show()
    
    
    
    
    