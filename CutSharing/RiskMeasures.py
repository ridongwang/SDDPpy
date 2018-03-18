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
    
    @abstractmethod
    def forward_pass_updates(self, *args, **kwargs):
        pass

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
    def forward_pass_updates(self, *args, **kwargs):
        pass    

class DistRobustDuality(AbstracRiskMeasure):
    INF_NORM = 'inf_norm'
    L1_NORM = 'L1_norm'
    L2_NORM = 'L2_norm'
    def __init__(self, dro_solver, dro_solver_params):
        super().__init__()
        self.cutting_planes_approx =  dro_solver_params['cutting_planes']
        self._dro_params = dro_solver_params
        
        if self.cutting_planes_approx == True:
            self.cuts_handler = None 
            self.cut_index = 0
        
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
            raise 'Risk measure does not support single cut'
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
            sp (StageProblem): Stage problem to modify.
            model (GRBModel): Model object associate to the stage.
            n_outcomes(int): Number of outcomes for the following stage.
            
        '''
        assert sp.multicut, 'This risk measure implementation is only compatible with multicut setting.'
        t = sp.stage
        if sp._last_stage:
            return
        
        
        set_type = self._dro_params['set_type']
        if set_type in [DistRobustDuality.L1_NORM, DistRobustDuality.L2_NORM]:
            r = self._dro_params['DUS_radius']
            q = self._dro_params['nominal_p']
        
            lambda_var =  model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='lambda[%i]' %(t))
            gamma_var =  model.addVar(lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='gamma[%i]' %(t))
            nu_var =  model.addVars(n_outcomes,  lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='nu[%i]' %(t))
            model.update()
            #Update objective function
            new_objective = sp.cx  + lambda_var + r*gamma_var - quicksum(q[i]*nu_var[i] for i in range(n_outcomes))
            model.setObjective(new_objective, GRB.MINIMIZE)
            model.addConstrs( (lambda_var - nu_var[i] - sp.oracle[i]>= 0 for i in range(n_outcomes)), 'dro_dual_ctr')
            if self.cutting_planes_approx == False:
                    
                if set_type == DistRobustDuality.L2_NORM:
                    model.params.QCPDual = 1
                    model.params.BarQCPConvTol = 1E-8
                    model.addConstr( quicksum(nu_var[i]*nu_var[i] for i in range(n_outcomes)) <= gamma_var*gamma_var    , 'norm_dro_ctr')
                elif set_type == DistRobustDuality.L1_NORM:
                    model.addConstrs((nu_var[i] + gamma_var>=0 for i in range(n_outcomes)) , 'norm_dro_ctr')
                    model.addConstrs((-nu_var[i] + gamma_var>=0 for i in range(n_outcomes)) , 'norm_dro_ctr')
            else: #Using cutting planes
                #Add max norm to bound
                model.addConstrs((nu_var[i] + gamma_var>=0 for i in range(n_outcomes)) , 'norm_dro_ctr')
                model.addConstrs((-nu_var[i] + gamma_var>=0 for i in range(n_outcomes)) , 'norm_dro_ctr')
                
                if set_type == DistRobustDuality.L2_NORM:
                    def g(x,order):
                        n = len(x)
                        gamma = x[-1]
                        if order == 0:
                            g_val = np.sqrt(sum(x[i]**2 for i in range(n-1))) - gamma
                            return g_val
                        elif order == 1:
                            norm_v = np.sqrt(sum(x[i]**2 for i in range(n-1)))
                            g_grad = [x[i]/norm_v for i in range(n)]
                            g_grad[-1] = -1
                            return g_grad
                        else:
                            raise 'Order of the function is either 0 or 1.'
                    var_names = [nu_var[v].VarName for v in nu_var]
                    var_names.append(gamma_var.VarName)
                    self.cuts_handler = DRO_CuttingPlanes([g] , var_names)
                    print(sp)    
                elif set_type == DistRobustDuality.L1_NORM:
                    self.cutting_planes_approx = False # There is nothing to approximate
                
            model.update()
        else:
            raise 'Set different from L1 and L2 norm are not supported'
            
    def forward_pass_updates(self, sp , **kwargs):
        '''
        Runs updates associated to the risk measure during the forward pass
        and after the model is optimized.
        '''
        
        tol = kwargs['fea_tol']
        if self.cutting_planes_approx == True and sp._last_stage == False:
            cph = self.cuts_handler
            var_vals = [sp.model.getVarByName(vn).X  for vn in cph.dual_set_var] 
            
            fea, ctr, vio = cph.check_dro_feasibility(var_vals, tol)
            if fea == False:
                model_vars = [sp.model.getVarByName(vn)  for vn in cph.dual_set_var] 
                cut_lhs = cph.refine_set(var_vals, model_vars,  ctr, vio)
                #print(sp.stage, '__>' , vio, cut_lhs)
                sp.model.addConstr(cut_lhs, GRB.LESS_EQUAL, 0, 'cut_dro_%i' %(self.cut_index))
                self.cut_index =  self.cut_index + 1
                return vio
        return 0

class DRO_CuttingPlanes():
    '''
    Implements helper function to do an outer approximation of the 
    of the uncertainty set in the DRO setting. In particular, it approximates 
    the dual a given uncertainty set
    
    Attrs:
        dual_dro_set (list of func): A list of functions that define the dual uncertatny set.    
            Each function has the signature:
                def g(vars_values, order):
                    vars_values: a point to evalute the function
                    order:  0 return the function evaluation
                            1 returns a subgradient in the same format as the input
            The set is assumed to be g_i(x,0) <= 0 \forall i
        dual_set_var (collection of str): A list of place holders for the variables that define the set
    '''
    
    def __init__(self, dual_dro_set , dual_dro_var):
        self.dual_dro_set = dual_dro_set
        self.dual_set_var = dual_dro_var
        
        
    def check_dro_feasibility(self, var_values, tolerance):
        '''
        Check feasibility of the set
        Args:
            var_values (collection of floats): Current value of the variables.
            tolerance (float): Feasibility tolerance.
        Return:
            status (bool): True is current point is feasible, False o.w.
            most_violated (int): index of the most volated inequality
            max_violation (float): violation of the constraint (value of the function)
        '''
        assert tolerance >=0, 'Feasibility tolerance must be positive.'
        assert type(var_values) == type(self.dual_set_var), 'Invalid format of the variables'
        most_violated = None
        max_violation = tolerance
        for (i,f) in enumerate(self.dual_dro_set):
            f_val = f(var_values, 0)
            if f_val>max_violation:
                max_violation = f_val
                most_violated = i
        
        if most_violated!=None:
            return False, most_violated, max_violation
        else:
            return True, None, 0
    
    def refine_set(self,var_values, vars,  g_func, g_val):
        '''
        Generate a valid inequality that refines the set
        Args:
            var_values (collection of floats): Current value of the variables.
            vars (collection of GRBVar): Variables of the model.
            g_func (func): Function associated to the most violated inequality.
            g_val (float): current value of the function
        Return:
            cut (LinExpr): A linear expression of the cut. 
        '''
        
        sub_gradient = self.dual_dro_set[g_func](var_values, 1)
        
        cut_intercept = 0 #  For 2 norm the intercept is zero! g_val - quicksum(sub_gradient[i]*var_values[i] for i in range(len(vars)))
        cut_lhs = quicksum(sub_gradient[i]*vars[i] for i in range(len(vars))) + cut_intercept
        return cut_lhs
        






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
        #p = srv.p
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

    def forward_pass_updates(self, *args, **kwargs):
        pass    


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
        #m.params.FeasibilityTol = 1E-9
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
        
        
        
        
        
        
        
