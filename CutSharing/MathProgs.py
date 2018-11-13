'''
Created on Nov 17, 2017

@author: dduque
'''
from gurobipy import GRB, Model, tupledict
from CutSharing.CutManagament import Cut
from CutSharing.CutManagament import CutPool
import CutSharing 
from CutSharing import *
from time import time
import scipy.sparse as sp
from CutSharing.RiskMeasures import Expectation

class StageProblem():
    '''
    Class to represent a general stage problem
    
    Attributes:
        model (GRB Model): Model of the stage
        cuts (CutPool): Object storing cuts
    '''


    def __init__(self, stage, model_builder, next_stage_rnd_vector, last_stage=False, risk_measure = Expectation(), multicut = False):
        '''
        Constructor
        
        Args:
            stage(int): Stage of the model
            model_builder(func): A function that build the stage model (GRBModel)
            next_stage_rnd_vector (StageRandomVector): Random vector that characterizes realizations in the next stage.
            last_stage (bool): Boolean flag indicating if this stage problem is the last one.
            risk_measure (AbstractRiskMeasure): A risk measure object to form the cuts.
            multicut (bool): Boolean flag indicating weather to use multicut or single cut.
        '''
        self.stage = stage
        self._last_stage = last_stage 
        self.model_stats = MathProgStats()
        model, in_states, out_states, rhs_vars = model_builder(stage)
        
        self.model = model
        self.risk_measure = risk_measure
        self.cut_pool = CutPool(stage)
        self.multicut = multicut
        
        self.states_map = {}  
        self.in_state = [x for x in in_states]
        self.out_state = [x for x in out_states]
        self.out_state_var  = {x:model.getVarByName(x) for x in out_states}
        self.gen_states_map(self.in_state)
        self.rhs_vars = [x for x in rhs_vars]
        self.rhs_vars_var  = {x:model.getVarByName(x) for x in rhs_vars}
        
        # Optimizer parameters
        self.model.params.OutputFlag = 0
        self.model.params.Threads = CutSharing.options['grb_threads']
        self.model.params.Method = 1
        #self.model.params.PreDual = 1
        #self.model.params.Presolve = 0
        #self.model.params.DualReductions = 0
        #self.model.params.FeasibilityTol = 1E-9
        
       
        # Add oracle var(s) and include it in the objective
        self.cx = self.model.getObjective() #Get objective before adding oracle variable
        if last_stage == False:
            num_outcomes = next_stage_rnd_vector.outcomes_dim 
            if multicut == False: #Gen single cut variable
                self.oracle = self.model.addVars(1,lb=-1E8, vtype = GRB.CONTINUOUS, name = 'oracle[%i]' %(stage))
            else:
                if num_outcomes == 0:
                    raise 'Multicut algorithm requires to define the number of outcomes in advance.'
                self.oracle = self.model.addVars(num_outcomes, lb=-1E8, vtype = GRB.CONTINUOUS, name = 'oracle[%i]' %(stage))
            self.model.update()
            risk_measure.modify_stage_problem(self, self.model, next_stage_rnd_vector)
              
                
        
        #Construct dictionaries of (constraints,variables) key where duals are needed
        self.ctrsForDuals = set() #Set of names
        self.ctrsForDualsRef = {} #Dict containing the reference to the constraints
        self.ctrInStateMatrix = {}
        for vname in self.in_state:
            var = self.model.getVarByName(vname)
            col =  self.model.getCol(var)
            for j in range(0,col.size()):
                ctr = col.getConstr(j)
                ctr_coeff = col.getCoeff(j)
                self.ctrsForDuals.add(ctr.ConstrName)
                self.ctrsForDualsRef[ctr.ConstrName]=ctr
                self.ctrInStateMatrix[ctr.ConstrName, vname] = -ctr_coeff
                
        self.ctrRHSvName = {}
        for vname in self.rhs_vars:
            var = self.model.getVarByName(vname)
            col =  self.model.getCol(var)
            #assert col.size() == 1, "RHS noise is not well defined"
            for c_index in range(col.size()):
                ctr = col.getConstr(c_index)
                assert ctr.ConstrName not in self.ctrRHSvName, "Duplicated RHS noise variable in the a single constraint."
                self.ctrRHSvName[ctr.ConstrName] = vname
                self.ctrsForDuals.add(ctr.ConstrName)
                self.ctrsForDualsRef[ctr.ConstrName]=ctr
        
           
                
     
    def solve(self, in_state_vals=None, random_realization=None, forwardpass = False, random_container=None, sample_path = None, num_cuts = 0):   
        '''
        Solves a stage problem given the state variables of the previous stage
        Args:
            in_state_vals (dict of real): dictionary containing the values of the state variables
                in the previous stage. They are referred as in_state for this model.
            random_realization (dic of real): dictionary containing the realization values to be solved.
                the key is the name of the variable that models the random value as a placeholder.
        
        Output:
            output (dict of objects): a dictionary with output information
                'status': Optimization model status
                'duals': dictionary of the dual variables associated to previous 
                        stage cuts.
                'out_state': Value of the state variables at the end of the sate
                        (input for the next stage).
                'cputime': float time solving the problem
                'lptipe': float value of the time solving the lp only
                'cutupdatetime' float value of the time updating the cuts.
        ''' 
        tnow = time()  
        cutupdatetime = 0  
        setuptime = 0;
        if self.stage > 0:
            #if forwardpass:
            setuptime = time()
            assert len(in_state_vals)==len(self.in_state), "In state vector has different cardinality than expected"
            for in_s_name in in_state_vals:
                in_s_val = in_state_vals[in_s_name]
                s_lb = self.out_state_var[in_s_name].lb
                s_ub = self.out_state_var[in_s_name].ub
                if in_s_val < s_lb:
                    in_s_val = s_lb
                elif in_s_val > s_ub: 
                    in_s_val = s_ub
                self.states_map[in_s_name].lb = in_s_val
                self.states_map[in_s_name].ub = in_s_val
            
            assert len(random_realization)==len(self.rhs_vars), "In random vector has different cardinality than expected %i" %(len(random_realization)-len(self.rhs_vars))
            for rr in random_realization:
                self.rhs_vars_var[rr].lb = random_realization[rr]
                self.rhs_vars_var[rr].ub = random_realization[rr]   
            self.model.update()
            setuptime = time()  - setuptime 
            
            cutupdatetime = time()   
            self.update_cut_pool(random_container, random_realization)
            cutupdatetime = time()  - cutupdatetime  
            
        
        #Solve LP
        lp_time = time()
        
        self.model.optimize()
      
        lp_time = time() - lp_time
        
        data_mgt_time = time()
        output = {}
        status = gurobiStatusCodeToStr(self.model.status)
        output['status'] = status
        if status ==  SP_OPTIMAL or status == SP_UNKNOWN:
            output['objval'] = self.model.objVal
                
            if forwardpass == True:
                resolves = 0
                resolve, violation = self.risk_measure.forward_pass_updates(self, fea_tol = 1E-2)
                while resolve and resolves<2 :
                    #print('Pass %i: resolving %i for violation %f' %(num_cuts, self.stage,violation))
                    lp_time_resolve = time()
                    self.model.optimize()
                    lp_time+= time() - lp_time_resolve
                    resolve, violation = self.risk_measure.forward_pass_updates(self, fea_tol = 1E-2)
                    resolves+=1
                    
                output['risk_measure_info'] = violation
                #if resolve == True and self.stage<=10000 and num_cuts > 50:
                    
                #    return self.solve(in_state_vals, random_realization, forwardpass, random_container, sample_path, num_cuts)
                #if resolve == True and self.stage<=3:
                #    print('Pass %i: resolving %i for violation %f' %(num_cuts, self.stage,violation))
                output['out_state'] = {vname:self.out_state_var[vname].X for vname in self.out_state_var}
            else:
                output['duals'] = {cname:self.ctrsForDualsRef[cname].Pi for cname in self.ctrsForDuals}
                output['cut_duals'] = {cut.name:cut.ctrRef.Pi for cut in self.cut_pool if cut.is_active}
                output['dual_obj_rhs_noice'] = sum(self.rhs_vars_var[self.ctrRHSvName[ctr_name]].UB*output['duals'][ctr_name] for ctr_name in self.ctrRHSvName)
            self.model_stats.add_simplex_iter_entr(self.model.IterCount)
        else:
            print('Not optimal sub')
            #if self.stage == 0:
                #self.print_theta()
            
        output['lptime'] = lp_time
        output['cutupdatetime'] = cutupdatetime
        output['setuptime']  = setuptime
        data_mgt_time = time() - data_mgt_time
        output['datamanagement'] = data_mgt_time
        output['cputime'] = time() - tnow
        
        return output
    

    
    def print_stage_res_summary(self):
        strout = ''
        #=======================================================================
        # for v in self.model.getVars():
        #     if 'generation' in v.varname or 'inflow[1' in v.varname or 'reservoir_level[1' in v.varname:
        #         strout = strout + '%20s:%10.3f;' %(v.varname, v.X)
        #=======================================================================
        strout = '======== Stage %s , HydroGen: %10.2f,  Thermal Gen %10.2f =========\n' %(self.stage,self.model.getVarByName('generation').X, self.model.getVarByName('thermal_gen').X)
        for r in range(10):
            v=self.model.getVarByName('inflow[%i,1]' %r)
            strout = strout + '+%10.3f; ' %(v.X)
        strout += '\n'
        for r in range(10):
            v=self.model.getVarByName('pour[%i]' %r)
            strout = strout + '+%10.3f; ' %(v.X)
        strout += '\n'
        for r in range(10):
            v=self.model.getVarByName('outflow[%i]' %r)
            strout = strout + '-%10.3f; ' %(v.X)
        strout += '\n'
        for r in range(10):
            v=self.model.getVarByName('spill[%i]' %r)
            strout = strout + '-%10.3f; ' %(v.X)
        strout += '\n'
        for r in range(10):
            v=self.model.getVarByName('reservoir_level[%i]' %r)
            strout = strout + '=%10.3f; ' %(v.X)
        print(strout)
    
    def print_theta(self):
        tvals = np.array([self.oracle[v].X for v in self.oracle])
        args_t_vals = np.argsort(tvals)
        print(args_t_vals)
            
        
    def printPostSolutionInformation(self):
        print('------------------------------------')
        print('Model in stage %i: obj>> %f' %(self.stage, self.model.ObjVal))
        for c in self.model.getConstrs():
            print('%20s %10.3f %10.3f'  %(c.ConstrName, c.RHS, c.PI))
        print('%20s %10s %10s %10s %10s %10s' %('name', 'lb' , 'ub' , 'obj', 'x', 'RC'))
        for v in self.model.getVars():
            print('%20s %10.3f %10.6e %10.3f %10.3f %10.3f'   %(v.varname, v.LB, v.UB, v.obj, v.X, v.RC))
        print('------------------------------------\n')
       
    #===========================================================================
    # def get_in_state_var(self, out_state):
    #     '''
    #     Returns the corresponding in state variable name
    #     given a out state variable name of the previous
    #     stage.
    #     '''
    #     sindex = out_state.index('[')
    #     newkey = out_state[:sindex]+'0'+out_state[sindex:]
    #     return newkey
    #===========================================================================
    def gen_states_map(self, in_states):
        for in_state in in_states:
            if '[' in in_state:
                sindex = in_state.index('[')
                #self.model.getVarByName(
                self.states_map[in_state] = self.model.getVarByName(in_state[:sindex-1]+in_state[sindex:])
                self.states_map[in_state[:sindex-1]+in_state[sindex:]] = self.model.getVarByName(in_state)
            else:
                sindex = in_state.index('0')
                self.states_map[in_state] = self.model.getVarByName(in_state[:sindex]+in_state[sindex+1:])
                self.states_map[in_state[:sindex]+in_state[sindex+1:]] = self.model.getVarByName(in_state)
                
        
    def get_out_state_var(self, in_state):
        '''
        Returns the corresponding out state variable name
        given a in state variable name of the next
        stage.
        '''
        return self.states_map[in_state].VarName
        #=======================================================================
        # sindex = in_state.index('[')
        # newkey = in_state[:sindex-1]+in_state[sindex:]
        # return newkey
        #=======================================================================

    
    def update_cut_pool(self, random_container, current_outcome):
        if self.stage > 0 and len(self.cut_pool)>0 and self.cut_pool.needs_update():
            srv = random_container[self.stage]
            omega_stage_abs_order = np.zeros((len(current_outcome),1))
            for rhs_c in current_outcome:
                omega_stage_abs_order[srv.vector_order[rhs_c]] = current_outcome[rhs_c]
            
            for cut in self.cut_pool:
                cut.adjust_intercept( omega_stage_abs_order )
            
            
    def createStageCut(self, cut_id, sp_next, rnd_vector_next, outputs_next, sample_path_forward_states , sample_path):
        '''
        Creates a cut for this stage problem
        
        Args:
            cut_id (int): Numeric id of the cut (corresponds to the number of backward passes).
            sp_next (StageProblem): stage problem object of the next stage
                TODO: This might not be enougth for lags>1
            rnd_vector_next (StageRandomVector): random vector containing the random variables 
                of the next stage.
            outputs_next (list of dict): list of outputs for every outcome in the next stage. 
                Each element of the list is a dictionary with the same structure as the output
                of ::func::StageProblem.solve.
            sample_path_forward_states(dict): dictionary with the values of the states of this
                stage that were computed in the forward pass associated with the current sample
                path.
            sample_path (list of dict): current sample path. Each element of the
                list corresponds to the realizations of a different stage.
        
        '''
        pi_bar = {}     #Expected duals of the transition function Ax = b + Bx0
        srv = rnd_vector_next
        soo = outputs_next
        spfs = sample_path_forward_states
        pi_bar, cut_gradiend_coeffs = self.risk_measure.compute_cut_gradient(self, sp_next, srv, soo, spfs , cut_id)
        cut_intercepts = self.risk_measure.compute_cut_intercept(self, sp_next, srv, soo, spfs, cut_id)
        
        stagewise_ind  = srv.is_independent
        
        if stagewise_ind:
            new_cuts = [Cut(self, cut_gradiend_coeffs[i],cut_intercepts[i], cut_id, outcome = i) for (i,grad) in enumerate(cut_gradiend_coeffs)]
            self.cut_pool.addCuts(self.model, new_cuts)
            #===================================================================
            # for (i,grad) in enumerate(cut_gradiend_coeffs):
            #     new_cut = Cut(self, cut_gradiend_coeffs[i],cut_intercepts[i], cut_id, outcome = i)
            #     self.cut_pool.addCut(self.model, new_cut)
            #===================================================================
        else:
            if self.multicut == True:
                raise "Multicut is not yet implemented for the dependent case"
            #ctrRHSvName
            alpha_bar = {}  #Expected duals of the cuts
            ab_D = np.zeros((1, len(pi_bar))) #Computation alpha_bar_{t+1}*D_{t+1}
            for cut in sp_next.cut_pool:
                alpha_bar[cut.name] = sum(srv.p[i]*soo[i]['cut_duals'][cut.name] for (i,o) in enumerate(srv.outcomes))
                ab_D =  ab_D + alpha_bar[cut.name]*cut.dep_rhs_vector
                
            omega_stage_abs_order = np.zeros((len(sample_path[self.stage]),1))
            pi_bar_abs_order = np.zeros((1, len(pi_bar)))
            for cpi in pi_bar:
                rhs_c = sp_next.ctrRHSvName[cpi]
                pi_bar_abs_order[0,srv.vector_order[rhs_c]] = pi_bar[cpi]
                omega_stage_abs_order[srv.vector_order[rhs_c]] = sample_path[self.stage][rhs_c]
            dep_rhs_vector = (sp.csr_matrix(pi_bar_abs_order + ab_D)).dot(srv.autoreg_matrices[-1]).toarray() #TODO: generalize for mor lags!!
            dep_rhs = dep_rhs_vector.dot(omega_stage_abs_order).squeeze()
            ind_rhs = cut_intercepts - dep_rhs
            new_cut = Cut(self, cut_gradiend_coeffs, cut_intercepts, cut_id, 
                          stagewise_ind=False,
                          ind_rhs=ind_rhs,
                          dep_rhs_vector=dep_rhs_vector)
            self.cut_pool.addCut(self.model, new_cut)
    
    def get_stage_objective_value(self):
        '''
        Returns the stage cost value
        '''
        try:
            return self.cx.getValue()
        except:
            #TODO: log for a warning
            print("ERROR IN get_stage_objective_value function")
            return 0
            
            
    
    def __repr__(self):
        return "SP(%i): #cuts:%i" %(self.stage,len(self.cut_pool.pool))

class MathProgStats:
    '''
    An object to keep track of some stats of the coresponding math program.
    
    Stats tracked:
        Mean number of iterations in simplex.
        Std of iteration in simplex.
    
    '''
    def __init__(self):
        self._simplex_iter_sum = 0.0
        self._simplex_iter_sumSq = 0.0
        self._simplex_iter_entries = 0
    
    def add_simplex_iter_entr(self, x):
        self._simplex_iter_entries+=1
        self._simplex_iter_sum += x
        self._simplex_iter_sumSq += (x**2)
    
    def get_mean(self):
        return self._simplex_iter_sum/self._simplex_iter_entries
    
    def get_sd(self):
        n =self._simplex_iter_entries
        numerator = self._simplex_iter_sumSq - (self._simplex_iter_sum**2)/n
        return np.sqrt(numerator/(n-1))
        