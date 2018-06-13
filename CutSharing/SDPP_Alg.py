'''
Created on Nov 17, 2017

@author: dduque
'''
import logging
import time

import numpy as np
import CutSharing as cs
from CutSharing.MathProgs import StageProblem, not_optimal_sp

from CutSharing.RandomManager import alg_rnd_gen, in_sample_gen, out_sample_gen, reset_all_rnd_gen

from OutputAnalysis.SimulationAnalysis import SimResult
from CutSharing.RandomnessHandler import ScenarioTree
from CutSharing.RiskMeasures import Expectation

 
sddp_log = cs.logger

alg_options = cs.alg_options()
iteration_log = ''
class SDDP(object):
    '''
    Implementation of Stochastic Dual Dynamic Programming algorithm.
    '''
    
    def __init__(self, T, model_builder, random_builder, risk_measure = Expectation, **risk_measure_params):
        '''
        Constructor
        '''
        self.stats = Stats()
        self.stage_problems = []
        self.random_container = random_builder()
        self.createStageProblems(T, model_builder, risk_measure, **risk_measure_params)
        
        self.instance = {'risk_measure':risk_measure, 'risk_measure_params':risk_measure_params}
        
        
        self.lb = None
        self.ub = float('inf')
        self.upper_bounds = [] 
        self.pass_iteration = 0
        self.num_cuts = 0
        
        #Attribute to keep track of the maximum violation when using cutting planes. 
        self.cutting_plane_max_vio = None

        
    def createStageProblems(self, T, model_builder, risk_measure, **risk_measure_params):
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
            sp_risk_measure = risk_measure(**risk_measure_params)
            next_stage_rnd_vector = self.random_container[i+1] if i<T-1 else None
            sp = StageProblem(i,model_builder, next_stage_rnd_vector,  i==T-1, risk_measure= sp_risk_measure, multicut = alg_options['multicut'])
            self.stage_problems.append(sp)

    def forwardpass(self, sample_path, simulation = False):
        '''
        Runs a forward pass given a sample path. If no sample pass is given,
        a dynamic version of the forward pass method is invoked. 
        '''
        if sample_path==None or len(sample_path)==0: #No sample path was given
            return self.dynamic_forwardpass(sample_path, simulation)
            
        fp_out_states = []
        fp_ub_value = 0
        for (i,sp) in enumerate(self.stage_problems):
            in_state = fp_out_states[-1] if i>0 else  None
            sp_output = sp.solve(in_state_vals = in_state, 
                                 random_realization= sample_path[i], 
                                 forwardpass = True, 
                                 random_container=self.random_container, 
                                 sample_path = sample_path,
                                 num_cuts = self.num_cuts)
            
            if sp_output['status'] != cs.SP_OPTIMAL:
                self.debrief_infeasible_sub(sample_path, i, sp_output, sp)
            
            fp_out_states.append(sp_output['out_state'])
            
            '''
            IO and stats updates
            '''
            if simulation and alg_options['outputlevel']>=3:
                    sp.print_stage_res_summary()
            if i == 0:       
                self.lb = sp_output['objval']
            fp_ub_value += sp.get_stage_objective_value()
            self.stats.updateStats(cs.FORWARD_PASS, lp_time=sp_output['lptime'], 
                                                         cut_update_time=sp_output['cutupdatetime'],
                                                         model_update_time=sp_output['setuptime'],
                                                         data_out_time= sp_output['datamanagement'],
                                                         num_lp_ctrs=sp.model.num_constrs,
                                                         iteration=self.pass_iteration)
            if sp_output['risk_measure_info']!=None and i==0:
                self.cutting_plane_max_vio =  sp_output['risk_measure_info']
                
        self.upper_bounds.append(fp_ub_value)
        if simulation and alg_options['outputlevel']>=3:
            print('---------------------------')
        return fp_out_states

    def dynamic_forwardpass(self, sample_path, simulation = False):
        '''
        Runs a forward pass given a sample path
        '''
        fp_out_states = []
        fp_ub_value = 0
        for (i,sp) in enumerate(self.stage_problems):
            in_state = fp_out_states[-1] if i>0 else  None
            self.random_container.getStageSample(i, sample_path, alg_rnd_gen)
            sp_output = sp.solve(in_state_vals = in_state, 
                                 random_realization= sample_path[i], 
                                 forwardpass = True, 
                                 random_container=self.random_container, 
                                 sample_path = sample_path,
                                 num_cuts = self.num_cuts)
            
            if sp_output['status'] != cs.SP_OPTIMAL:
                self.debrief_infeasible_sub(sample_path, i, sp_output, sp)
            
            fp_out_states.append(sp_output['out_state'])
            sp.risk_measure.forward_prob_update(i,self.random_container)   
            '''
            IO and stats updates
            '''
            if simulation and alg_options['outputlevel']>=3:
                    sp.print_stage_res_summary()
            if i == 0:       
                self.lb = sp_output['objval']
            fp_ub_value += sp.get_stage_objective_value()
            self.stats.updateStats(cs.FORWARD_PASS, lp_time=sp_output['lptime'], 
                                                         cut_update_time=sp_output['cutupdatetime'],
                                                         model_update_time=sp_output['setuptime'],
                                                         data_out_time= sp_output['datamanagement'],
                                                         num_lp_ctrs=sp.model.num_constrs,
                                                         iteration=self.pass_iteration)
            if sp_output['risk_measure_info']!=None and i==0:
                self.cutting_plane_max_vio =  sp_output['risk_measure_info']
                
        self.upper_bounds.append(fp_ub_value)
        if simulation and alg_options['outputlevel']>=3:
            print('---------------------------')
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
                                                         data_out_time= sp_output['datamanagement'],
                                                         num_lp_ctrs=sp.model.num_constrs,
                                                         iteration=self.pass_iteration)
                
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
        
    def init_out(self, instance_name):
        sddp_log.info(instance_name)
        sddp_log.info('T: %i' %(len(self.stage_problems)))
        sddp_log.info('Number of passes: %i:' %(alg_options['max_iter']))
        sddp_log.info('Sample paths per pass: %i:' %(alg_options['n_sample_paths']))
        sddp_log.info('Multicut: %s' %(str(alg_options['multicut'])))
        
        if alg_options['outputlevel']>=2:
            sddp_log.info('==============================================================================================')
            sddp_log.info('%4s %15s %15s %15s %12s %15s'
                  %('Pass', 'LB', 'iUB', 'iHW', 'Wall time', 'Other'))
            sddp_log.info('==============================================================================================')
        
    
    def iteration_update(self,fp_time, bp_time, force_print = False, last_iter = False):
        if (alg_options['outputlevel']>=2 and (self.pass_iteration % alg_options['lines_freq'] == 0 or last_iter==True) and  force_print==False):
            elapsed_time = time.time() - self.ini_time
            additional_msg = '' 
            if self.cutting_plane_max_vio !=None:
                additional_msg = '%15.5e' %(self.cutting_plane_max_vio)
            sddp_log.info('%4i %15.5e %15.5e %15.5e %15.2f %15s' %(self.pass_iteration, self.lb, self.ub, self.ub_hw, elapsed_time, additional_msg))
        if  force_print:
            sddp_log.info('==============================================================================================')
            sddp_log.info('%4s %15s %15s %15s %12s %12s %12s'
                          %('Pass', 'LB', 'iUB','iHW', 'F time', 'B time', 'Wall time'))
            elapsed_time = time.time() - self.ini_time
            sddp_log.info('%4s %15.5e %15.5e %15.5e %12.2f %12.2f %12.2f' 
                  %("Sim%i" %(alg_options['sim_iter']) , self.lb, self.ub, self.ub_hw, fp_time, bp_time, elapsed_time))
            sddp_log.info('==============================================================================================')
        
    def termination(self):
        if self.pass_iteration >= alg_options['max_iter']:
            return True
        if self.pass_iteration > 0:
            if self.lb >= self.ub - self.ub_hw - alg_options['opt_tol']:
                return True             
        return False
    
    def run(self, pre_sample_paths = None, ev = False, instance_name = 'Default', dynamic_sampling = False):
        reset_all_rnd_gen()
        lbs = []
        self.ini_time = time.time()
        self.init_out(instance_name)
        fp_time = 0
        bp_time = 0
        termination = False
        while termination ==False:
            
            f_timer = time.time()
            sample_paths = []
            fp_outputs = [] 
            
            for i in range(0,alg_options['n_sample_paths']):
                s_path = None
                if pre_sample_paths == None and dynamic_sampling == False:
                    s_path, _ = self.random_container.getSamplePath(alg_rnd_gen, ev = ev)
                elif pre_sample_paths != None:
                    s_path = pre_sample_paths.pop()
                else:
                    s_path = list()
                
                output_fp = self.forwardpass(sample_path = s_path)
                sample_paths.append(s_path)
                fp_outputs.append(output_fp)
            fp_time = time.time()-f_timer
            self.stats.updateStats(cs.FORWARD_PASS, total_time = fp_time)
            
            '''
            Compute statistical upper bounds
            '''
            if self.pass_iteration % 10 == 0:
                self.compute_statistical_bound(alg_options['in_sample_ub'])
                if self.pass_iteration>0:
                    self.compute_upper_bound_math_prog(30)
            
            '''
            Stopping criteria
            '''
            lbs.append(self.lb)
            termination =  self.termination()
            self.iteration_update(fp_time, bp_time, last_iter = termination)
            if termination:
                break
            
            b_timer = time.time()
            for i in range(0,alg_options['n_sample_paths']):
                s_path = sample_paths[i]
                output_fp = fp_outputs[i]
                self.backwardpass(forward_out_states = output_fp, sample_path=s_path, ev=ev)
            bp_time = time.time()-b_timer
            self.stats.updateStats(cs.BACKWARD_PASS, total_time = bp_time)
        
            self.pass_iteration+=1
            
            
        #self.stats.print_report(instance_name, self.stage_problems)
        return lbs
    
    
    def simulate_policy(self, n_samples, out_of_sample_random_container):
        np.random.seed(1111)
        self.upper_bounds = []
        for i in range(0,n_samples):
            s_path, _ =  out_of_sample_random_container.getSamplePath(out_sample_gen)
            #===================================================================
            # if i<=2:
            #         print(s_path[3]['innovations[6]'],' ', s_path[9]['innovations[5]'])
            #===================================================================
            if alg_options['outputlevel']>=3:
                sddp_log.debug('Simulation %i:' %(i))
                sddp_log.debug(s_path)
            
            output_fp = self.forwardpass(sample_path = s_path, simulation=True)
        self.ub = np.mean(self.upper_bounds)
        self.ub_hw = 2*np.std(self.upper_bounds)/np.sqrt(len(self.upper_bounds))
        sr = SimResult(self.instance, self.upper_bounds.copy())
        if alg_options['outputlevel']>=1:
            self.iteration_update(0, 0, force_print = True)
        return sr
    
    def compute_statistical_bound1(self, n_samples):
        self.upper_bounds = []  # rest bound
        for i in range(0,n_samples):
            s_path, _ =  self.random_container.getSamplePath(in_sample_gen)
            if alg_options['outputlevel']>=3:
                sddp_log.debug('Simulation %i:' %(i))
                sddp_log.debug(s_path)
            output_fp = self.forwardpass(sample_path = s_path, simulation=True)
        
        self.ub = np.mean(self.upper_bounds)
        self.ub_hw = 2*np.std(self.upper_bounds)/np.sqrt(len(self.upper_bounds))
            
    def compute_statistical_bound(self, n_samples):
        '''
        Computes an statistical upper bound. For DRO risk measures,
        probabilities in the tree change for every sample path to account for
        the worse case expectation. Both the upper bound and halfwidth are stored
        as an attribute of the class.
        
        Args:
            n_samples(int): Number of sample paths to compute the upper bound
        
        NOTES: For DRO risk measure, the worst case probabilities are achieved
                asymptotically with the number of iterations of SDDP. 
        '''
        self.upper_bounds = []  # rest bound
        for k in range(0,n_samples):
            fp_out_states = []
            fp_ub_value = 0
            sample_path = [] #partial sample path
            for (t,sp) in enumerate(self.stage_problems):
                self.random_container.getStageSample(t, sample_path, in_sample_gen)
                in_state = fp_out_states[-1] if t>0 else  None
                sp_output = sp.solve(in_state_vals = in_state, 
                                     random_realization= sample_path[t], 
                                     forwardpass = True, 
                                     random_container=self.random_container, 
                                     sample_path = sample_path,
                                     num_cuts = self.num_cuts)
                if sp_output['status'] != cs.SP_OPTIMAL:
                    self.debrief_infeasible_sub(sample_path, t,sp_output,sp)
                
                fp_out_states.append(sp_output['out_state'])
                fp_ub_value += sp.get_stage_objective_value()
                
                sp.risk_measure.forward_prob_update(t,self.random_container)   
            self.upper_bounds.append(fp_ub_value)     
        self.ub = np.mean(self.upper_bounds)
        self.ub_hw = 2*np.std(self.upper_bounds)/np.sqrt(len(self.upper_bounds))
        self.random_container.reset_to_nominal_dist()
        
    def compute_upper_bound_math_prog(self, n_samples):
        '''
        Computes an upper bound assigning probabilities to sample paths' solutions.
        First, all sample paths are realized and solved. These sample paths induce a 
        subtree of the full scenario tree which is used to assign probabilities to each 
        sample path. An auxiliary linear program is build to found the maximum (worst-case)
        expectation of the realized samples, introducing constraints to enforce that the 
        distribution at every stage is within the uncertainty set.  
        
        Args:
            n_samples(int): Number of sample paths to compute the upper bound
        '''
        sub_tree = ScenarioTree(self.random_container)
        for k in range(0,n_samples):
            fp_out_states = []
            fp_ub_value = 0
            sample_path = [] #partial sample path
            sample_path_outcomes = [] #partial sample path
            for (t,sp) in enumerate(self.stage_problems):
                _, stage_outcome=  self.random_container.getStageSample(t, sample_path, in_sample_gen)
                sample_path_outcomes.append(stage_outcome)
                in_state = fp_out_states[-1] if t>0 else  None
                sp_output = sp.solve(in_state_vals = in_state, 
                                     random_realization= sample_path[t], 
                                     forwardpass = True, 
                                     random_container=self.random_container, 
                                     sample_path = sample_path,
                                     num_cuts = self.num_cuts)
                if sp_output['status'] != cs.SP_OPTIMAL:
                    self.debrief_infeasible_sub(sample_path, t,sp_output,sp)
                
                fp_out_states.append(sp_output['out_state'])
                fp_ub_value += sp.get_stage_objective_value()
            
            sub_tree.add_sample_path(k, sample_path, sample_path_outcomes, fp_ub_value)
 
        ub  = sub_tree.compute_subtree_upper_bound([sp.risk_measure for sp in self.stage_problems])
        print(ub)
        
        
        
        
        
    def debrief_infeasible_sub(self, sample_path, t,sp_output,sp): 
        for ss in sample_path: 
            print(ss)
            print('GRB STATUS %i' %(sp.model.status))
            sp.model.write('%s%i_.lp' %(sp_output['status'], t))
            sp.model.computeIIS()
            sp.model.write("model.ilp")
            raise not_optimal_sp('A stage %i problem was not optimal' %(t))
        
class Stats:
    '''
    Class to keep track of performance statistics.
    '''
    def __init__(self):
        self.lp_time = {cs.FORWARD_PASS:0.0, cs.BACKWARD_PASS:0.0}
        self.cut_update_time = {cs.FORWARD_PASS:0.0, cs.BACKWARD_PASS:0.0}
        self.cut_gen_time = {cs.FORWARD_PASS:0.0, cs.BACKWARD_PASS:0.0}
        self.pass_time = {cs.FORWARD_PASS:0.0, cs.BACKWARD_PASS:0.0}
        self.model_update_time = {cs.FORWARD_PASS:0.0, cs.BACKWARD_PASS:0.0}
        self.lp_counter = {cs.FORWARD_PASS:0, cs.BACKWARD_PASS:0}
        self.data_out = {cs.FORWARD_PASS:0.0, cs.BACKWARD_PASS:0.0}
        self.lp_times = []
        
    def updateStats(self, passType=cs.FORWARD_PASS,
                     lp_time = 0.0, 
                     cut_update_time = 0.0, 
                     cut_gen_time=0.0, 
                     model_update_time=0.0, 
                     data_out_time = 0.0,
                     total_time = 0.0,
                     num_lp_ctrs = 0,
                     iteration = 0):
        self.lp_time[passType] += lp_time
        self.pass_time[passType] += total_time
        self.cut_update_time[passType] += cut_update_time
        self.cut_gen_time[passType] += cut_gen_time
        self.model_update_time[passType] += model_update_time
        self.data_out[passType] += data_out_time
        if lp_time > cs.ZERO_TOL:
            self.lp_counter[passType] += 1
            self.lp_times.append( (lp_time , num_lp_ctrs , iteration) )
    
    def print_lp_data(self, instance_name,stage_problems):
        sddp_log.info('Simplex Iterations Stats')
        sddp_log.info('%10s %12s %12s %12s' %('Stage', 'Mean # iter', 'SD # iter', '# Entries'))
        for sp in stage_problems:
            stage = sp.stage
            mean_iter = sp.model_stats.get_mean()
            sd_iter = sp.model_stats.get_sd()
            n_entries = sp.model_stats._simplex_iter_entries
            sddp_log.info('%10i %12.2f %12.2f %12.2f' %(stage, mean_iter, sd_iter, n_entries))
            
            
        
        import csv
        file_name = '%s/%s.csv' %(cs.cwd,instance_name)
        with open(file_name, 'w') as myfile:
            fieldnames = ['pass', 'num_ctr' , 'lp_time']
            writer = csv.DictWriter(myfile, fieldnames=fieldnames)
            writer.writeheader()
            for x in self.lp_times:    
                writer.writerow({'lp_time':"%f" %(x[0]), 'num_ctr':"%i" %(x[1]), 'pass':"%i" %(x[2])})    
    
    
    def print_report(self,instance_name , stage_problems):  
        #self.print_lp_data(instance_name, stage_problems)
        sddp_log.info('Time profiling')
        sddp_log.info('%15s %12s %12s %12s %12s %12s %12s %12s %12s' %('Pass', '# LPs', 'setup', 'simplex', 'output', 'cut update' ,'cut gen', 'other','total'))
        
        ''' FORWARD PASS '''
        f_mu = self.model_update_time[cs.FORWARD_PASS]
        f_lp = self.lp_time[cs.FORWARD_PASS]
        f_cu = self.cut_update_time[cs.FORWARD_PASS]
        f_cg = self.cut_gen_time[cs.FORWARD_PASS]
        f_do = self.data_out[cs.FORWARD_PASS]
        f_tot = self.pass_time[cs.FORWARD_PASS]
        f_other = f_tot - (f_mu+ f_lp+ f_cu+f_cg + f_do)
        f_lp_count = self.lp_counter[cs.FORWARD_PASS]
        sddp_log.info('%15s %12i %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f' \
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
        sddp_log.info('%15s %12i %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f' \
              % ('Backward',b_lp_count, b_mu, b_lp, b_do, b_cu, b_cg, b_other, b_tot))
