'''
Created on Nov 13, 2017

@author: dduque
'''
import numpy as np

class ScenarioSet():
    
    def __init__(self, scenarioDim, sce_list ):
        self.dim = scenarioDim
        self.scenarios = sce_list
        self.expectation = None
     
    def checkProbConsisntency(self,):
        prob = sum(self.scenarios[k].p for k in self.scenarios)
        if np.abs(prob  - 1.0) < 1E-8:
            return True
        else:
            return False
           
    def computeExpectation(self):
        self.expectation = (sum(self.scenarios[k].p*self.scenarios[k].scenarioVector for k in self.scenarios))
        return self.expectation
        
    def computeCovarianve(self):
        weightedVectors = np.zeros((self.dim, len(self.scenarios)))
        i = 0 
        probs = [sc.p for sc in self.scenarios.values()] 
        for sc in self.scenarios.values():
            weightedVectors[:,i] = sc.scenarioVector
            i = i+1
            
            
        
        return np.cov(weightedVectors, aweights = probs )
       
        
class Scenario():
    '''
    Class to represent a scenario of a random vector
   '''
    def __init__(self, nid, scenarioValues, probability):
        '''
        Attributes:
            nid (int): unique id of the scenario
            scenarioVector (List of numbers): values that define a scenario
        '''
        assert isinstance(scenarioValues, np.ndarray)
        self.id = nid
        self.scenarioVector = scenarioValues
        self.p = probability
    
    def computeDistance(self, distFunc, other):  
        '''
        Args:
         distFunc (:func:): Function to compute the distance between two scenarios.
         other (Scenario): Another scenario to measure the distance wrt to it.
        '''
        value = distFunc(self,other)
        return value

    def findClosest(self, distFunc, reducedSet):
        '''
            Finds the closest scenario in the reduce set
            
            Return:
                Scenarion in the armin_k distFunc(self,reduceSet[k])
        '''
        minValue = float("inf")
        minScenario = None
        for k in reducedSet:
            d = self.computeDistance(distFunc, reducedSet[k])
            if d<minValue:
                minValue = d
                minScenario = reducedSet[k]
        return minScenario
    def __repr__(self):
        return '%i' %(self.id) 