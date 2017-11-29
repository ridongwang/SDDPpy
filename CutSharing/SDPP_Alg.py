'''
Created on Nov 17, 2017

@author: dduque
'''
from CutSharing.MathProgs import StageProblem,not_optimal_sp
import CutSharing
import numpy as np

alg_options = CutSharing.alg_options()

class SDDP(object):
    '''
    classdocs
    '''

    def __init__(self, T, model_builder, random_builder):
        '''
        Constructor
        '''
        self.stage_problems = []
        self.createStageProblems(T, model_builder)
        self.random_container = random_builder()
        
        self.random_container.preprocess_randomness()
        
        self.lb = None
        self.ub = float('inf')
        self.upper_bounds = [] 
        self.pass_iteration = 0
        
        
    def createStageProblems(self, T, model_builder):
        '''
        Creates all subproblems given a builder.
        Args:
            T (int): Number of stages.
            model_builder (::func::) a function the return a math model.
            in_states (list of str): list with the names of the variables that 
                represent the previous state.
            out_states (list of str): list with the names of the variables that 
                represent the next state.
        '''
        for i in range(T):
            sp = StageProblem(i,model_builder, i==T-1)
            self.stage_problems.append(sp)
   
    def termination(self):
        if self.pass_iteration > 1:
            ub_bar = np.mean(self.upper_bounds)
            ub_sd = np.std(self.upper_bounds)
            if np.abs(self.lb  - ub_bar - 2*ub_sd/np.sqrt(len(self.upper_bounds)))<CutSharing.SDDP_OPT_TOL:
                return True
            self.ub = ub_bar
                
        if self.pass_iteration == alg_options['max_iter']:
            return True
        return False

    
    def forwardpass(self, sample_path):
        '''
        Runs a forward pass given a sample path
        '''
        fp_out_states = []
        fp_ub_value = 0
        for (i,sp) in enumerate(self.stage_problems):
            in_state = fp_out_states[-1] if i>0 else  None
            sp_output = sp.solve(in_state_vals = in_state, 
                                 random_realization= sample_path[i], 
                                 forwardpass = True, 
                                 random_container=self.random_container, 
                                 sample_path = sample_path)
            
            #assert sp_output['status'] == CutSharing.SP_OPTIMAL, sp_output['status']
            if sp_output['status'] != CutSharing.SP_OPTIMAL:
                for ss in sample_path: print(ss);
                raise not_optimal_sp('A stage problem was not optimal')
            fp_out_states.append(sp_output['out_state'])
            if i == 0:
                #sp.printPostSolutionInformation()
                self.lb = sp_output['objval']
            fp_ub_value += sp.get_stage_objective_value()
        self.upper_bounds.append(fp_ub_value)
        return fp_out_states
    
    def backwardpass(self, forward_out_states = None,  sample_path = None):
        '''
        Runs a backward pass given a sample path and the forward pass decision
        associated to the sample path.
        
        Args:
            forward_out_states (list of dict): List of the out states for each stage.
            The out states are the solutions of decision variables in the problem and
            are stored as a dictionary where the key is the variable name.
            
            sample_path (list of outcomes): List of the sample path being solved
        '''
        T = len(self.stage_problems)
        
        for t in range(T-1, 0, -1):
            outputs_per_outcome = []
            sp = self.stage_problems[t]
            stage_rnd_vector = self.random_container[t]
            omega_t = stage_rnd_vector.getOutcomes(sample_path)
            for outcome in omega_t:
                sp_output = sp.solve(in_state_vals=None, 
                                     random_realization=outcome, 
                                     forwardpass=False, 
                                     random_container=self.random_container, 
                                     sample_path = sample_path)
                
                #sp.printPostSolutionInformation()
                outputs_per_outcome.append(sp_output)
            self.createStageCut(t-1, stage_rnd_vector, outputs_per_outcome, omega_t, sample_path)
            del(outputs_per_outcome)
        return None
        
    def createStageCut(self, stage, stage_rnd_vector, outputs_per_outcome, omega, sample_path):
        '''
        Calls the cut routine in the stage given as a parameter.
        Args:
            stage (int): Stage where the cut will be added
            stage_rnd_vector (StageRandomVector): object with the stochastic 
                representation of the next stage.
            outputs_per_outcome (list of dict): a list with the outputs of the
                of the next stage (one per scenario of the current pass).
            omega (list of dict): a list of scenarios of the current sample 
                path for the next stage. Each scenario is represented as a 
                dictionary where the key is the random element name and the
                value is the the realization in each scenario.
            sample_path (list of dict): current sample path. Each element of the
                list corresponds to the realizations of a different stage.
        '''
        if stage<0:
            return
        sp_t = self.stage_problems[stage]
        sp_t1 = self.stage_problems[stage+1]
        sp_t.createStageCut(self.pass_iteration, sp_t1, stage_rnd_vector, outputs_per_outcome,omega, sample_path )  
        
        
    def iteration_update(self):
        
        if alg_options['outputlevel']>=2:
            print('%3i %12.5e %12.5e' %(self.pass_iteration, self.lb, self.ub))
        self.pass_iteration+=1
    
    def run(self):
        while self.termination()==False:
            s_path = self.random_container.getSamplePath()
            output_fp = self.forwardpass(sample_path = s_path)
            output_bp = self.backwardpass(forward_out_states = output_fp, sample_path=s_path)
            self.iteration_update()
        
