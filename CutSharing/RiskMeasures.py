'''
Created on Jan 11, 2018

@author: dduque
'''
from abc import ABC, abstractmethod
from gurobipy import *
import CutSharing as cs
import numpy as np
import logging
sddp_log = cs.logger

class AbstracRiskMeasure(ABC):
    '''
    classdocs
    -Define what a risk measure should have
    -Relate this to the random handler clasess. 
    
    '''
    
    
    @abstractmethod
    def __init__(self):
        '''
        Constructor
        '''
        self._static_dist = False
        self._current_cut_gradient = None
    
    @abstractmethod
    def compute_cut_gradient(self):
        raise('Method not implemented in abstract class')
    @abstractmethod
    def compute_cut_intercept(self):
        raise('Method not implemented in abstract class')
    

class Expectation(AbstracRiskMeasure):
    
    def __init__(self):
        super().__init__()
        
    def compute_cut_gradient(self, sp, sp_next, srv, soo, spfs):
        '''
        Args:
            
        '''
        pi_bar = {}
        for ctr in sp_next.ctrsForDuals:
            pi_bar[ctr] = sum(srv.p[i]*soo[i]['duals'][ctr] for (i,o) in enumerate(srv.outcomes))
            
        cut_gradiend_coeff = {vo:0 for vo in sp.out_state}
        for (c,vi) in sp_next.ctrInStateMatrix:
            vo = sp_next.get_out_state_var(vi)
            cut_gradiend_coeff[vo] += pi_bar[c]*sp_next.ctrInStateMatrix[c,vi]
        
        self._current_cut_gradient = cut_gradiend_coeff
        return pi_bar, cut_gradiend_coeff
    
    def compute_cut_intercept(self, sp, sp_next, srv, soo, spfs):
        cut_gradiend_coeff = self._current_cut_gradient
        cut_intercept = sum(srv.p[i]*soo[i]['objval'] for (i,o) in enumerate(srv.outcomes)) - sum(spfs[vn]*cut_gradiend_coeff[vn] for vn in sp.out_state) 
        return cut_intercept






'''
Distributionaly robust classes
'''
class DistRobust(AbstracRiskMeasure):
    
    def __init__(self, dro_solver, dro_solver_params):
        '''
        inner_solver (DRO_solver): a class that computes the worst
            distribution to update the risk measure probabilities.
        '''
        super().__init__()
        self.inner_solver = dro_solver(**dro_solver_params)
        self._wors_case_dist = None
        self._static_dist = True
       
        
    def compute_cut_gradient(self, sp, sp_next, srv, soo, spfv):
        '''
        Args:
            
        '''
        zs = np.array([soo[i]['objval'] for i in range(len(srv.outcomes))])
        p = self.inner_solver.compute_worst_case_distribution(zs)
        self._wors_case_dist = p
        pi_bar = {}
        for ctr in sp_next.ctrsForDuals:
            pi_bar[ctr] = sum(p[i]*soo[i]['duals'][ctr] for (i,o) in enumerate(srv.outcomes))
            
        cut_gradiend_coeff = {vo:0 for vo in sp.out_state}
        for (c,vi) in sp_next.ctrInStateMatrix:
            vo = sp_next.get_out_state_var(vi)
            cut_gradiend_coeff[vo] += pi_bar[c]*sp_next.ctrInStateMatrix[c,vi]
        self._current_cut_gradient =cut_gradiend_coeff
        return pi_bar, cut_gradiend_coeff
    
    def compute_cut_intercept(self, sp, sp_next, srv, soo, spfs):
        cut_gradiend_coeff = self._current_cut_gradient
        p = self._wors_case_dist
        cut_intercept = sum(p[i]*soo[i]['objval'] for (i,o) in enumerate(srv.outcomes)) - sum(spfs[vn]*cut_gradiend_coeff[vn] for vn in sp.out_state) 
        return cut_intercept
        
        
        
class DistRobusInnerSolver(ABC):
    @abstractmethod
    def __init__(self):
        pass
    
    @abstractmethod
    def compute_worst_case_distribution(self):
        pass

class PhilpottInnerDROSolver(DistRobusInnerSolver):
    def __init__(self, nominal_p, DUS_radius):
        self.nominal_p = nominal_p
        self.uncertanty_radius = DUS_radius
        self._one_time_warning = True
        
        self.model = self.build_model(nominal_p, DUS_radius)
    
    def build_model(self, q, DUS_radius):
        m = Model('DRO_solver')
        m.params.OutputFlag = 0 
        p = m.addVars(len(q) , lb = 0, ub = 1, obj=1, vtype=GRB.CONTINUOUS, name='p')
        m.update()
        m.addConstr(p.sum(), GRB.EQUAL, 1, 'prob')
        m.addConstr(quicksum((p[i]- q[i])*(p[i]- q[i]) for i in p), GRB.LESS_EQUAL, DUS_radius**2, 'DUS')
        m.setObjective(m.getObjective(), GRB.MAXIMIZE)
        m.update()
        return m
    def compute_worst_case_distribution(self, outcomes_objs):
        '''
        Compute the worst cas probability distribution for a particular stage
        given the objective function value of the descendent nodes (as many as outcomes)
        '''
        vars = self.model.getVars()
        for (i,z) in enumerate(outcomes_objs):
            vars[i].obj = z
        self.model.update()
        self.model.optimize()
        
        
        if self._one_time_warning:
            self._one_time_warning = False
            sddp_log.warning('DRO inner problem not implemented!')
        assert len(outcomes_objs) == len(self.nominal_p)
        new_p = np.array([vars[i].X for i in range(len(vars))])
        return new_p
        
        
        
        
        
        
        