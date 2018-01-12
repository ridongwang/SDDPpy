'''
Created on Jan 11, 2018

@author: dduque
'''

class AbstracRiskMeasure(object):
    '''
    classdocs
    '''
    
    -Define what a risk measure should have
    -Relate this to the random handler clasess. 
    

    def __init__(self):
        '''
        Constructor
        '''
    
        self._static_dist = False
    

class Expectation(AbstracRiskMeasure):
    
    pass

class DistRobust(AbstracRiskMeasure):
    
    def __init__(self, inner_solver):
        '''
        inner_solver (DRO_solver): a class that computes the worst
            distribution to update the risk measure probabilities.
        '''
        super().__init__()
        self._static_dist = True
        for i in ran