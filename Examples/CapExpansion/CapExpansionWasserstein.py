'''
Created on Apr 3, 2019

@author: dduque
'''
import numpy as np
from gurobipy import GRB, Model
from CutSharing.RandomnessHandler import RandomContainer, StageRandomVector
from Examples.CapExpansion.CapExpDataGen import Datahandler
from CutSharing.SDDP_Alg import SDDP
from Utils.file_savers import write_object_results 
from CutSharing.SDDP_utils import report_stats
import CutSharing
from CutSharing.RiskMeasures import DiscreteWassersteinInnerSolver, DistRobust,\
    DistRobustWasserstein
import pickle
from OutputAnalysis.SimulationAnalysis import plot_sim_results
import os

def random_bulder_factory(data_cap_exp):
    def random_builder():
        rc = RandomContainer()
        for t in range(2):
            rv_t = StageRandomVector(t)
            rc.append(rv_t)
            for j in data_cap_exp.J:
                if t>0:
                    rv_t.addRandomElememnt('d[%i]' %(j), data_cap_exp.d[j,:])
                else:
                    rv_t.addRandomElememnt('d[%i]' %(j), [0.0])
        return rc

    return random_builder;
def stage_prob_factory(data):
    def stage_prob_builder(t):
        model =  Model('CapExp%i' %t)   #Create a model
        x = model.addVars(data.I,vtype=GRB.CONTINUOUS,lb=0,ub=data.b,obj=0,name='x');
        x0 = model.addVars(data.I,vtype=GRB.CONTINUOUS,lb=0,ub=data.b,obj=0,name='x0');
        d  = model.addVars(data.J,vtype=GRB.CONTINUOUS,lb=0,ub=0,obj=0,name='d');
        model.update()
        if t == 0:
            model.addConstr(sum(x[i] for i in data.I), GRB.LESS_EQUAL, data.b, 'globalCap')
            model.update()
        else:
            y = model.addVars(data.I,
                                    data.J,
                                    vtype=GRB.CONTINUOUS, 
                                    lb=0,
                                    obj=data.c,
                                    name='y');
            s = model.addVars(data.J,
                                    vtype=GRB.CONTINUOUS, 
                                    lb=0,
                                    obj=data.rho,
                                    name='s');
            model.update()
            '''Capacity constraints''' 
            model.addConstrs((sum(y[i,j] for j in data.J)==x0[i] for i in data.I), 'CapCtr');
            '''Recourse constraint, subcontracting '''
            model.addConstrs((sum(y[i,j] for i in data.I)+s[j]  >= d[j] for j in data.J), 'RecCtr');
            model.update()
    
        in_states = [var.VarName for var in x0.values()]
        out_states = [var.VarName for var in x.values()]
        rhs_vars = [var.VarName for var in d.values()]
        return model, in_states, out_states, rhs_vars 
    return stage_prob_builder

if __name__ == '__main__':
    CutSharing.options['max_iter'] = 50
    CutSharing.options['multicut'] = True
    CutSharing.options['lines_freq'] = 1
    T=2 # number of stages
    ''' 
    Test problem
    Parameters
    '''
    m = 20  # number of facilities
    n = 100  # number of customers
    B = 7.5 * n  # budget
    rho = 10  # shortfall penalty
    
    scen = None # number of scenarios
    for scen in [200]:#[5,10,30,50,100]:
        '''
            Data generation
        '''
        data = Datahandler(n, m, scen, B, rho, seed=0) 
        data_oos = Datahandler(n, m, 1000, B, rho, seed=123456, dh=data) 
        model_builder = stage_prob_factory(data)
        random_builder = random_bulder_factory(data)
        dro_trials = [r for r in range(10,100,10)]
        dro_trials.extend([r for r in range(100,200,10)])
        #dro_trials.extend([r for r in range(200,1000,50)])
        #dro_trials.extend([r for r in range(1000,2001,200)])
        dro_trials.extend([0.5,1,2,5,200,300,400,500,1000,5000])
        dro_trials.sort()
        #dro_trials = [10000000]
        print(dro_trials)
      
        if False: #solve
            for dro_r in dro_trials:
                instance_name = 'DW_CapExp_N%i_m%i_n%i_r%.3f' %(scen,m,n,dro_r)
                #algo_sp = SDDP(T, model_builder, random_builder)
                algo_dro = SDDP(T, model_builder, random_builder, risk_measure = DistRobustWasserstein , norm = 1 , radius = dro_r, data_random_container=random_builder())
                #algo_dro = SDDP(T, model_builder, random_builder, risk_measure = DistRobust, dro_inner_solver = DiscreteWassersteinInnerSolver,  radius = dro_r, norm=1, data_random_container=random_builder())
                algo = algo_dro
                lbs = algo.run(instance_name=instance_name)
                
                save_path = './Output/%s_LBS.pickle' %(instance_name)
                write_object_results(save_path, (algo.instance, lbs))      
                
                random_builder_oos  = random_bulder_factory(data_oos)
                sim_result = algo.simulate_policy(random_builder_oos())
                save_path = './Output/%s_OOS.pickle' %(instance_name)
                write_object_results(save_path, sim_result)
                report_stats(sim_result.sims_ub) 
    
            # Run SP   
            instance_name = 'SP_CapExp_N%i_m%i_n%i' %(scen,m,n)
            algo_sp = SDDP(T, model_builder, random_builder)
            algo = algo_sp
            lbs = algo.run(instance_name=instance_name)
            
            save_path = './Output/%s_LBS.pickle' %(instance_name)
            write_object_results(save_path, (algo.instance, lbs))      
            
            random_builder_oos  = random_bulder_factory(data_oos)
            sim_result = algo.simulate_policy(random_builder_oos())
            save_path = './Output/%s_OOS.pickle' %(instance_name)
            write_object_results(save_path, sim_result)
            report_stats(sim_result.sims_ub)   
        
        
        if True: # Plot
            ptf='/Users/dduque/Dropbox/WORKSPACE/SDDP/Examples/CapExpansion/Output/'
            experiment_files = os.listdir(ptf)
            print(experiment_files)
            sim_results = []
            for dro_r in dro_trials:
                instance_name = 'DW_CapExp_N%i_m%i_n%i_r%.3f' %(scen,m,n,dro_r)
                new_sim = pickle.load(open('%s%s_OOS.pickle' %(ptf,instance_name), 'rb'))
                sim_results.append(new_sim)
        
            instance_name = 'SP_CapExp_N%i_m%i_n%i' %(scen,m,n)
            sp_sim = pickle.load(open('%s%s_OOS.pickle' %(ptf,instance_name), 'rb'))
            
            #Sort experiments
            if 'radius' in sim_results[0].instance['risk_measure_params']:
                sim_results.sort(key= lambda x:x.instance['risk_measure_params']['radius'])
            elif 'dro_solver_params' in  sim_results[0].instance['risk_measure_params']:
                if 'radius' in sim_results[0].instance['risk_measure_params']['dro_solver_params']:
                    sim_results.sort(key= lambda x:x.instance['risk_measure_params']['dro_solver_params']['radius'])
                elif 'DUS_radius' in sim_results[0].instance['risk_measure_params']['dro_solver_params']:
                    sim_results.sort(key= lambda x:x.instance['risk_measure_params']['dro_solver_params']['DUS_radius'])
                else:
                    print(sim_results[0].instance)
                    raise 'Unknown dro params'
            else:
                print(sim_results[0].instance)
                raise 'Unknown dro params'
            instance_name = 'CapExp_N%i_m%i_n%i' %(scen,m,n)
            plot_path = '%s%s.pdf' %(ptf,instance_name)
            N = scen
            plot_sim_results(sp_sim, sim_results, plot_path, N, excel_file=False)
