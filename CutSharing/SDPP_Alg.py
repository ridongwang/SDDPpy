'''
Created on Nov 17, 2017

@author: dduque
'''
from CutSharing.MathProgs import StageProblem,not_optimal_sp
import CutSharing as cs
import numpy as np
import time
from astropy import samp

alg_options = cs.alg_options()

class SDDP(object):
    '''
    classdocs
    '''

    def __init__(self, T, model_builder, random_builder):
        '''
        Constructor
        '''
        self.stats = Stats()
        self.stage_problems = []
        self.createStageProblems(T, model_builder)
        self.random_container = random_builder()
        
        self.random_container.preprocess_randomness()
        
        self.lb = None
        self.ub = float('inf')
        self.upper_bounds = [] 
        self.pass_iteration = 0
        self.num_cuts = 0
        
        
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

    def forwardpass(self, sample_path , simulation = False):
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
            
            
            
            assert sp_output['status'] == cs.SP_OPTIMAL, sp_output['status']
            if sp_output['status'] != cs.SP_OPTIMAL:
                for ss in sample_path: print(ss);
                print('GRB STATUS %i' %(sp.model.status))
                sp.model.write('%s%i_.lp' %(sp_output['status'], i))
                sp.model.computeIIS()
                sp.model.write("model.ilp")
                raise not_optimal_sp('A stage %i problem was not optimal' %(i))
            fp_out_states.append(sp_output['out_state'])
            if simulation and alg_options['outputlevel']>=3:
                    sp.print_stage_res_summary()
            if i == 0:       
                self.lb = sp_output['objval']
            fp_ub_value += sp.get_stage_objective_value()
            self.stats.updateStats(cs.FORWARD_PASS, lp_time=sp_output['lptime'], 
                                                         cut_update_time=sp_output['cutupdatetime'],
                                                         model_update_time=sp_output['setuptime'],
                                                         data_out_time= sp_output['datamanagement'])
        self.upper_bounds.append(fp_ub_value)
        
        return fp_out_states
    
    def backwardpass(self, forward_out_states = None,  sample_path = None, ev = False):
        '''
        Runs a backward pass given a sample path and the forward pass decision
        associated to the sample path.
        
        Args:
            forward_out_states (list of dict): List of the out states for each stage.
            The out states are the solutions of decision variables in the problem and
            are stored as a dictionary where the key is the variable name.
            
            sample_path (list of outcomes): List of the sample path being solved
            
            ev (bool): If expected value policy is being computed (Default is False).
        '''
        T = len(self.stage_problems)
        
        for t in range(T-1, 0, -1):
            outputs_per_outcome = []
            sp = self.stage_problems[t]
            stage_rnd_vector = self.random_container[t]
            omega_t = stage_rnd_vector.getOutcomes(sample_path, ev)
            for outcome in omega_t:
                sp_output = sp.solve(in_state_vals=forward_out_states[t-1], 
                                     random_realization=outcome, 
                                     forwardpass=False, 
                                     random_container=self.random_container, 
                                     sample_path = sample_path)
                
                #sp.printPostSolutionInformation()
                outputs_per_outcome.append(sp_output)
                self.stats.updateStats(cs.BACKWARD_PASS, lp_time=sp_output['lptime'], 
                                                         cut_update_time=sp_output['cutupdatetime'],
                                                         model_update_time=sp_output['setuptime'],
                                                         data_out_time= sp_output['datamanagement'])
                
            cut_creation_time = self.createStageCut(t-1, stage_rnd_vector, outputs_per_outcome, forward_out_states[t-1], sample_path)
            self.stats.updateStats(cs.BACKWARD_PASS, cut_gen_time=cut_creation_time)
            #del(outputs_per_outcome)
        self.num_cuts +=1   
        
        
    def createStageCut(self, stage, stage_rnd_vector, outputs_per_outcome, forward_out_states, sample_path):
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
            forward_out_states 
            sample_path (list of dict): current sample path. Each element of the
                list corresponds to the realizations of a different stage.
            
        '''
        cut_creation_time = time.time()
        if stage<0:
            return
        sp_t = self.stage_problems[stage]
        sp_t1 = self.stage_problems[stage+1]
        sp_t.createStageCut(self.num_cuts, sp_t1, stage_rnd_vector, outputs_per_outcome, forward_out_states, sample_path )  
        return time.time()-cut_creation_time
        
    def printInit(self):
        print('T: %i' %(len(self.stage_problems)))
        print('Number of passes: %i:' %(alg_options['max_iter']))
        print('Sample paths per pass: %i:' %(alg_options['n_sample_paths']))
        print('==========================================================================================')
        print('%3s %15s %15s %15s %12s %12s %12s'
              %('Pass', 'LB', 'UB^','UB_hw', 'F time', 'B time', 'Wall time'))
        print('==========================================================================================')
            
    def iteration_update(self,fp_time, bp_time, force_print = False):
        self.ub = np.mean(self.upper_bounds)
        self.ub_hw = 2*np.std(self.upper_bounds)/len(self.upper_bounds)
        if (alg_options['outputlevel']>=2 and self.pass_iteration % alg_options['lines_freq'] == 0 and force_print==False):
            elapsed_time = time.time() - self.ini_time
            print('%3i %15.5e %15.5e %15.5e %12.2f %12.2f %12.2f' 
                  %(self.pass_iteration, self.lb, self.ub, self.ub_hw,fp_time, bp_time, elapsed_time))
        if  force_print:
            elapsed_time = time.time() - self.ini_time
            print('==========================================================================================')
            print('%3s %15.5e %15.5e %15.5e %12.2f %12.2f %12.2f' 
                  %("Sim" , self.lb, self.ub, self.ub_hw,fp_time, bp_time, elapsed_time))
    def termination(self):
        if self.pass_iteration >= alg_options['max_iter']-1:
            return True
        if self.pass_iteration > 0:
            if (self.ub - self.ub_hw) < self.lb  < (self.ub +self.ub_hw):
                return False                
        return False
    
    def run(self, pre_sample_paths = None, ev = False):
        lbs = []
        self.ini_time = time.time()
        self.printInit()
        fp_time = 0
        bp_time = 0
        termination = False
        while termination ==False:
            
            f_timer = time.time()
            sample_paths = []
            fp_outputs = [] 
            self.upper_bounds = [] #rest upper bounds
            for i in range(0,alg_options['n_sample_paths']):
                s_path = None
                if pre_sample_paths == None:
                    s_path = self.random_container.getSamplePath( ev = ev)
                else:
                    s_path = pre_sample_paths.pop()
                output_fp = self.forwardpass(sample_path = s_path)
                sample_paths.append(s_path)
                fp_outputs.append(output_fp)
            fp_time = time.time()-f_timer
            self.stats.updateStats(cs.FORWARD_PASS, total_time = fp_time)
            
            '''
            Stopping criteria
            '''
            lbs.append(self.lb)
            self.iteration_update(fp_time, bp_time)
            if self.termination():
                termination = True
                break
            
            b_timer = time.time()
            for i in range(0,alg_options['n_sample_paths']):
                s_path = sample_paths[i]
                output_fp = fp_outputs[i]
                self.backwardpass(forward_out_states = output_fp, sample_path=s_path, ev=ev)
            bp_time = time.time()-b_timer
            self.stats.updateStats(cs.BACKWARD_PASS, total_time = bp_time)
        
            self.pass_iteration+=1
            
            
        self.stats.printReport()
        return(lbs)
    
    
    def simulate_policy(self, n_samples):
        self.upper_bounds = []
        for i in range(0,n_samples):
            s_path =  self.random_container.getSamplePath()
            if alg_options['outputlevel']>=3:
                print('Simulation %i:' %(i))
                print(s_path)
            
            output_fp = self.forwardpass(sample_path = s_path, simulation=True)
            
        self.iteration_update(0, 0, force_print = True)
        
class Stats:
    
    def __init__(self):
        self.lp_time = {cs.FORWARD_PASS:0.0, cs.BACKWARD_PASS:0.0}
        self.cut_update_time = {cs.FORWARD_PASS:0.0, cs.BACKWARD_PASS:0.0}
        self.cut_gen_time = {cs.FORWARD_PASS:0.0, cs.BACKWARD_PASS:0.0}
        self.pass_time = {cs.FORWARD_PASS:0.0, cs.BACKWARD_PASS:0.0}
        self.model_update_time = {cs.FORWARD_PASS:0.0, cs.BACKWARD_PASS:0.0}
        self.lp_counter = {cs.FORWARD_PASS:0, cs.BACKWARD_PASS:0}
        self.data_out = {cs.FORWARD_PASS:0.0, cs.BACKWARD_PASS:0.0}
        
        
    def updateStats(self, passType=cs.FORWARD_PASS,
                     lp_time = 0.0, 
                     cut_update_time = 0.0, 
                     cut_gen_time=0.0, 
                     model_update_time=0.0, 
                     data_out_time = 0.0,
                     total_time = 0.0):
        self.lp_time[passType] += lp_time
        self.pass_time[passType] += total_time
        self.cut_update_time[passType] += cut_update_time
        self.cut_gen_time[passType] += cut_gen_time
        self.model_update_time[passType] += model_update_time
        self.data_out[passType] += data_out_time
        if lp_time > cs.ZERO_TOL:
            self.lp_counter[passType] += 1
    
    def printReport(self):  
        print('Time profiling')
        print('%15s %12s %12s %12s %12s %12s %12s %12s %12s' %('Pass', '# LPs', 'setup', 'simplex', 'output', 'cut update' ,'cut gen', 'other','total'))
        
        ''' FORWARD PASS '''
        f_mu = self.model_update_time[cs.FORWARD_PASS]
        f_lp = self.lp_time[cs.FORWARD_PASS]
        f_cu = self.cut_update_time[cs.FORWARD_PASS]
        f_cg = self.cut_gen_time[cs.FORWARD_PASS]
        f_do = self.data_out[cs.FORWARD_PASS]
        f_tot = self.pass_time[cs.FORWARD_PASS]
        f_other = f_tot - (f_mu+ f_lp+ f_cu+f_cg + f_do)
        f_lp_count = self.lp_counter[cs.FORWARD_PASS]
        print('%15s %12i %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f' \
              % ('Forward', f_lp_count, f_mu, f_lp, f_do, f_cu, f_cg, f_other, f_tot))
        
        ''' BACKWARD PASS '''
        b_mu = self.model_update_time[cs.BACKWARD_PASS]
        b_lp = self.lp_time[cs.BACKWARD_PASS]
        b_cu = self.cut_update_time[cs.BACKWARD_PASS]
        b_cg = self.cut_gen_time[cs.BACKWARD_PASS]
        b_do = self.data_out[cs.BACKWARD_PASS]
        b_tot = self.pass_time[cs.BACKWARD_PASS]
        b_other = b_tot - (b_mu + b_lp + b_cu + b_cg + b_do)
        b_lp_count = self.lp_counter[cs.BACKWARD_PASS]
        print('%15s %12i %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f' \
              % ('Backward',b_lp_count, b_mu, b_lp, b_do, b_cu, b_cg, b_other, b_tot))
