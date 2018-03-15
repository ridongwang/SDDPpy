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
    Abstract representation of what a risk measure should do in SDD
    '''
    
    @abstractmethod
    def __init__(self):
        '''
        
        Attr:
            _current_cut_gradient (dict): local copy of the current
                cut gradient (used again when computing the.
        '''
        
        
        self._current_cut_gradient = None
    
    @abstractmethod
    def compute_cut_gradient(self):
        raise('Method not implemented in abstract class')
    @abstractmethod
    def compute_cut_intercept(self):
        raise('Method not implemented in abstract class')
    @abstractmethod
    def update_cut_intercept(self):
        raise('Method not implemented in abstract class')
    
    @abstractmethod
    def modify_stage_problem(self, *args):
        raise('Method not implemented in abstract class')

class Expectation(AbstracRiskMeasure):
    
    def __init__(self):
        super().__init__()
        
    def compute_cut_gradient(self, sp, sp_next, srv, soo, spfs):
        '''
        Computes expected dual variables for the single cut version
        and then the gradient.
        
        Args:
            sp (StageProblem): current subproblem where the cut will be added.
            sp_nex (StageProblem): subproblem for the next stage.
            srv (StageRandomVector): Random vector of the next stage
            soo (List of dict): A list of outputs of all the subproblems descendants.
            spfs (dict (str-float)): Values for the states of the current stage computed 
                in the forward pass.
        Return 
            pi_bar (list[dict]): Expected value of the duals. For the single cut algorithm
                the list contains just one element.
            cut_gradiend_coeff(list[dict]):
        '''
        multicut = sp.multicut
        pi_bar = None
        cut_gradiend_coeff = None
        
        if multicut == False: #single cut
            pi_bar = [{}]
            cut_gradiend_coeff = [{vo:0 for vo in sp.out_state}]
            for ctr in sp_next.ctrsForDuals:
                pi_bar[0][ctr] = sum(srv.p[i]*soo[i]['duals'][ctr] for (i,o) in enumerate(srv.outcomes))
            
            for (c,vi) in sp_next.ctrInStateMatrix:
                vo = sp_next.get_out_state_var(vi)
                cut_gradiend_coeff[0][vo] += pi_bar[0][c]*sp_next.ctrInStateMatrix[c,vi]
        else:           #Multicut
            pi_bar = [{} for _ in srv.outcomes]
            cut_gradiend_coeff = [{vo:0 for vo in sp.out_state} for _ in srv.outcomes]
            for ctr in sp_next.ctrsForDuals:
                for (i,o) in enumerate(srv.outcomes):
                    pi_bar[i][ctr] = srv.p[i]*soo[i]['duals'][ctr]
            
            for (c,vi) in sp_next.ctrInStateMatrix:
                vo = sp_next.get_out_state_var(vi)
                for (i,o) in enumerate(srv.outcomes):
                    cut_gradiend_coeff[i][vo] += pi_bar[i][c]*sp_next.ctrInStateMatrix[c,vi]
            
        self._current_cut_gradient = cut_gradiend_coeff
        return pi_bar, cut_gradiend_coeff
    
    def compute_cut_intercept(self, sp, sp_next, srv, soo, spfs):
        '''
        Computes cut intercept(s) 
        
        Args: (as in cut gradient method)
        Returns:
            cut_intercepts (list of int): List of intercep(s) of the cut(s).
        '''
        cut_gradiend_coeff = self._current_cut_gradient
        cut_intercepts = None
        if sp.multicut == False:
            cut_intercepts = [sum(srv.p[i]*soo[i]['objval'] for (i,o) in enumerate(srv.outcomes)) - sum(spfs[vn]*cut_gradiend_coeff[0][vn] for vn in sp.out_state)]
        else:
            cut_intercepts = [0  for _ in srv.outcomes]
            for (i,o) in enumerate(srv.outcomes):
                cut_intercepts[i] = srv.p[i]*soo[i]['objval'] - sum(spfs[vn]*cut_gradiend_coeff[i][vn] for vn in sp.out_state)
        return cut_intercepts

    
    def update_cut_intercept(self):
        pass
    def modify_stage_problem(self, sp,  model, n_outcomes):
        model.setObjective(sp.cx + sp.oracle.sum())
            

class DistRobustDuality(AbstracRiskMeasure):
    INF_NORM = 'inf_norm'
    L1_NORM = 'L1_norm'
    L2_NORM = 'L2_norm'
    def __init__(self,dro_solver, dro_solver_params):
        super().__init__()
        self._dro_params = dro_solver_params
    def compute_cut_gradient(self, sp, sp_next, srv, soo, spfs):
        '''
        Computes expected dual variables for the single cut version
        and then the gradient.
        
        Args:
            sp (StageProblem): current subproblem where the cut will be added.
            sp_nex (StageProblem): subproblem for the next stage.
            srv (StageRandomVector): Random vector of the next stage
            soo (List of dict): A list of outputs of all the subproblems descendants.
            spfs (dict (str-float)): Values for the states of the current stage computed 
                in the forward pass.
        Return 
            pi_bar (list[dict]): Expected value of the duals. For the single cut algorithm
                the list contains just one element.
            cut_gradiend_coeff(list[dict]):
        '''
        multicut = sp.multicut
        pi_bar = None
        cut_gradiend_coeff = None
        
        if multicut == False: #single cut
            raise 'Risk measure does not suppor single cut'
        else:           #Multicut
            pi_bar = [{} for _ in srv.outcomes]
            cut_gradiend_coeff = [{vo:0 for vo in sp.out_state} for _ in srv.outcomes]
            for ctr in sp_next.ctrsForDuals:
                for (i,o) in enumerate(srv.outcomes):
                    pi_bar[i][ctr] = soo[i]['duals'][ctr]
            
            for (c,vi) in sp_next.ctrInStateMatrix:
                vo = sp_next.get_out_state_var(vi)
                for (i,o) in enumerate(srv.outcomes):
                    cut_gradiend_coeff[i][vo] += pi_bar[i][c]*sp_next.ctrInStateMatrix[c,vi]
            
        self._current_cut_gradient = cut_gradiend_coeff
        return pi_bar, cut_gradiend_coeff
    
    def compute_cut_intercept(self, sp, sp_next, srv, soo, spfs):
        '''
        Computes cut intercept(s) 
        
        Args: (as in cut gradient method)
        Returns:
            cut_intercepts (list of int): List of intercep(s) of the cut(s).
        '''
        cut_gradiend_coeff = self._current_cut_gradient
        cut_intercepts = None
        if sp.multicut == False:
            cut_intercepts = [sum(soo[i]['objval'] for (i,o) in enumerate(srv.outcomes)) - sum(spfs[vn]*cut_gradiend_coeff[0][vn] for vn in sp.out_state)]
        else:
            cut_intercepts = [0  for _ in srv.outcomes]
            for (i,o) in enumerate(srv.outcomes):
                cut_intercepts[i] = soo[i]['objval'] - sum(spfs[vn]*cut_gradiend_coeff[i][vn] for vn in sp.out_state)
        return cut_intercepts
    
    def update_cut_intercept(self):
        pass
    
    def modify_stage_problem(self, sp,  model, n_outcomes):
        '''
        Modify the stage problem model to incorporate the DRO
        risk measure as dual variables of the inner problem.
        
        Args:
            sp (StageProblem): Stage problem to modify
            
        '''
        assert sp.multicut, 'This risk measure implementation is only compatible with multicut setting.'
        t = sp.stage
        if sp._last_stage:
            return
        
        r = self._dro_params['DUS_radius']
        q = self._dro_params['nominal_p']
        set_type = self._dro_params['set_type']
        lambda_var =  model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='lambda[%i]' %(t))
        gamma_var =  model.addVar(lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='gamma[%i]' %(t))
        nu_var =  model.addVars(n_outcomes,  lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='nu[%i]' %(t))
        model.update()
        #Update objective function
        new_objective = sp.cx  + lambda_var + r*gamma_var - quicksum(q[i]*nu_var[i] for i in range(n_outcomes))
        model.setObjective(new_objective, GRB.MINIMIZE)
        model.addConstrs( (lambda_var - nu_var[i] - sp.oracle[i]>= 0 for i in range(n_outcomes)), 'dro_dual_ctr')
        if set_type == DistRobustDuality.L2_NORM:
            model.params.QCPDual = 1
            model.params.BarQCPConvTol = 1E-8
            model.addConstr( quicksum(nu_var[i]*nu_var[i] for i in range(n_outcomes)) <= gamma_var*gamma_var    , 'norm_dro_ctr')
        elif set_type == DistRobustDuality.L1_NORM:
            model.addConstrs((nu_var[i] + gamma_var>=0 for i in range(n_outcomes)) , 'norm_dro_ctr')
            model.addConstrs((-nu_var[i] + gamma_var>=0 for i in range(n_outcomes)) , 'norm_dro_ctr')
        
        model.update()
        
        
        
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
        Computes expected dual variables for the single cut version
        and then the gradient.
        
        Args:
            sp (StageProblem): current subproblem where the cut will be added.
            sp_nex (StageProblem): subproblem for the next stage.
            srv (StageRandomVector): Random vector of the next stage
            soo (List of dict): A list of outputs of all the subproblems descendants.
            spfs (dict (str-float)): Values for the states of the current stage computed 
                in the forward pass.
        Return 
            pi_bar (list[dict]): Expected value of the duals. For the single cut algorithm
                the list contains just one element.
            cut_gradiend_coeff(list[dict]):
        '''
        zs = np.array([soo[i]['objval'] for i in range(len(srv.outcomes))])
        p = self.inner_solver.compute_worst_case_distribution(zs)
        self._wors_case_dist = p
        #=======================================================================
        # pi_bar = {}
        # for ctr in sp_next.ctrsForDuals:
        #     pi_bar[ctr] = sum(p[i]*soo[i]['duals'][ctr] for (i,o) in enumerate(srv.outcomes))
        #     
        # cut_gradiend_coeff = {vo:0 for vo in sp.out_state}
        # for (c,vi) in sp_next.ctrInStateMatrix:
        #     vo = sp_next.get_out_state_var(vi)
        #     cut_gradiend_coeff[vo] += pi_bar[c]*sp_next.ctrInStateMatrix[c,vi]
        # self._current_cut_gradient =cut_gradiend_coeff
        # return pi_bar, cut_gradiend_coeff
        #=======================================================================
        multicut = sp.multicut
        pi_bar = None
        cut_gradiend_coeff = None
        
        if multicut == False: #single cut
            pi_bar = [{}]
            cut_gradiend_coeff = [{vo:0 for vo in sp.out_state}]
            for ctr in sp_next.ctrsForDuals:
                pi_bar[0][ctr] = sum(p[i]*soo[i]['duals'][ctr] for (i,o) in enumerate(srv.outcomes))
            
            for (c,vi) in sp_next.ctrInStateMatrix:
                vo = sp_next.get_out_state_var(vi)
                cut_gradiend_coeff[0][vo] += pi_bar[0][c]*sp_next.ctrInStateMatrix[c,vi]
        else:           #Multicut
            pi_bar = [{} for _ in srv.outcomes]
            cut_gradiend_coeff = [{vo:0 for vo in sp.out_state} for _ in srv.outcomes]
            for ctr in sp_next.ctrsForDuals:
                for (i,o) in enumerate(srv.outcomes):
                    pi_bar[i][ctr] = p[i]*soo[i]['duals'][ctr]
            
            for (c,vi) in sp_next.ctrInStateMatrix:
                vo = sp_next.get_out_state_var(vi)
                for (i,o) in enumerate(srv.outcomes):
                    cut_gradiend_coeff[i][vo] += pi_bar[i][c]*sp_next.ctrInStateMatrix[c,vi]
            
        self._current_cut_gradient = cut_gradiend_coeff
        return pi_bar, cut_gradiend_coeff
    
    def compute_cut_intercept(self, sp, sp_next, srv, soo, spfs):
        '''
        Computes cut intercept(s) 
        
        Args: (as in cut gradient method)
        Returns:
            cut_intercepts (list of int): List of intercep(s) of the cut(s).
        '''
        #=======================================================================
        # cut_gradiend_coeff = self._current_cut_gradient
        # p = self._wors_case_dist
        # cut_intercept = sum(p[i]*soo[i]['objval'] for (i,o) in enumerate(srv.outcomes)) - sum(spfs[vn]*cut_gradiend_coeff[vn] for vn in sp.out_state) 
        # return cut_intercept
        #=======================================================================
        p = self._wors_case_dist
        cut_gradiend_coeff = self._current_cut_gradient
        cut_intercepts = None
        if sp.multicut == False:
            cut_intercepts = [sum(p[i]*soo[i]['objval'] for (i,o) in enumerate(srv.outcomes)) - sum(spfs[vn]*cut_gradiend_coeff[0][vn] for vn in sp.out_state)]
        else:
            cut_intercepts = [0  for _ in srv.outcomes]
            for (i,o) in enumerate(srv.outcomes):
                cut_intercepts[i] = p[i]*soo[i]['objval'] - sum(spfs[vn]*cut_gradiend_coeff[i][vn] for vn in sp.out_state)
        return cut_intercepts
        
        
        
    def update_cut_intercept(self):
        pass  
    def modify_stage_problem(self,  sp,  model, n_outcomes):
        model.setObjective(sp.cx + sp.oracle.sum())





class DistRobusInnerSolver(ABC):
    @abstractmethod
    def __init__(self):
        pass
    
    @abstractmethod
    def compute_worst_case_distribution(self):
        pass

class PhilpottInnerDROSolver(DistRobusInnerSolver):
    def __init__(self, nominal_p, DUS_radius, set_type):
        self.nominal_p = nominal_p
        self.uncertanty_radius = DUS_radius
        self._one_time_warning = True
        
        self.model, self.p_var = self.build_model(nominal_p, DUS_radius, set_type)
    
    def build_model(self, q, DUS_radius, set_type):
        m = Model('DRO_solver')
        m.params.OutputFlag = 0 
        p = m.addVars(len(q) , lb = 0, ub = 1, obj=1, vtype=GRB.CONTINUOUS, name='p')
        m.update()
        m.addConstr(p.sum(), GRB.EQUAL, 1, 'prob')
        if set_type == DistRobustDuality.L2_NORM:
            m.addConstr(quicksum((p[i]- q[i])*(p[i]- q[i]) for i in p), GRB.LESS_EQUAL, DUS_radius**2, 'DUS')
        elif set_type == DistRobustDuality.L1_NORM:
            d_plus = m.addVars(len(q) , lb = 0, ub = 1, obj=0, vtype=GRB.CONTINUOUS, name='dp')
            d_minus = m.addVars(len(q) , lb = 0, ub = 1, obj=0, vtype=GRB.CONTINUOUS, name='dm')
            m.addConstr(quicksum(d_plus[i] + d_minus[i] for i in p), GRB.LESS_EQUAL, DUS_radius, 'DUS')
            m.addConstrs((d_plus[i] - d_minus[i] == p[i]- q[i] for i in p), 'ABS_liner')
        m.setObjective(m.getObjective(), GRB.MAXIMIZE)
        m.update()
        return m, p
    def compute_worst_case_distribution(self, outcomes_objs):
        '''
        Compute the worst cas probability distribution for a particular stage
        given the objective function value of the descendent nodes (as many as outcomes)
        Args:
            outcomes_objs (ndarray): vector of objectives values of the next stage.
        Returns:
            new_p (ndarray): array with the worse case probabilities in DRO
        '''
        vars = self.model.getVars()
        for (i,z) in enumerate(outcomes_objs):
            vars[i].obj = z
        self.model.update()
        self.model.optimize()
        
        if self._one_time_warning:
            self._one_time_warning = False
            #sddp_log.warning('DRO inner problem not implemented!')
        assert len(outcomes_objs) == len(self.nominal_p)
        new_p = np.array([vars[i].X for i in range(len(vars))])
        return new_p
        
        
        
        
        
        
        