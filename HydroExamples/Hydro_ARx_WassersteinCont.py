'''
Created on Jan 3, 2018

@author: dduque
'''
import csv
import CutSharing
import logging
import numpy as np
from Utils.file_savers import write_object_results
np.set_printoptions(linewidth= 200, nanstr='nen')
from CutSharing.RandomnessHandler import RandomContainer, StageRandomVector, AR1_depedency
from CutSharing.SDPP_Alg import SDDP
from CutSharing import logger as sddp_log

from Utils.argv_parser import sys,parse_args
from gurobipy import *
from InstanceGen.ReservoirChainGen import read_instance, HydroRndInstance
from HydroExamples import *
from CutSharing.RiskMeasures import DistRobust, PhilpottInnerDROSolver, DistRobustDuality,\
    InnerDROSolverX2, DistRobustWasserstein, mod_chi2, DistRobustWassersteinCont
from OutputAnalysis.SimulationAnalysis import plot_sim_results,\
    plot_metrics_comparison, plot_lbs
from CutSharing.RandomManager import experiment_desing_gen,\
    reset_experiment_desing_gen

'''
Global variables to store instance data
'''
T = None
nr = None
lag = None
dro_radius = None
Rmatrix = None
RHSnoise = None
initial_inflow = None 
valley_chain = None
valley_chain_oos = None
prices = None
Water_Penalty = 10000

def random_builder():
    rc = RandomContainer()
    rndVectors = []
    for t in range(0,T):
        rv_t = StageRandomVector(t)
        rc.append(rv_t)
        for (i,r) in enumerate(valley_chain):
            if t>0:
                re = rv_t.addRandomElememnt('innovations[%i]' %(i), r.inflows)
            else:
                re = rv_t.addRandomElememnt('innovations[%i]' %(i), [0.0])
            rndVectors.append(rv_t)
    rc.preprocess_randomness()
    return rc

def random_builder_out_of_sample(vally_chain_sample):
    '''
    Generates a random container for out-of-sample performance.
    '''
    rc = RandomContainer()
    rndVectors = []
    for t in range(0,T):
        rv_t = StageRandomVector(t)
        rc.append(rv_t)
        for (i,r) in enumerate(vally_chain_sample):
            if t>0:
                re = rv_t.addRandomElememnt('innovations[%i]' %(i), r.inflows)
            else:
                re = rv_t.addRandomElememnt('innovations[%i]' %(i), [0.0])
            rndVectors.append(rv_t)
    rc.preprocess_randomness()
    return rc



def model_builder(stage):
    '''
    Builds a particular instance of a multistage problem
    '''
    import gurobipy as gb
    m = Model('Hydrovalley')
    
    '''
    State variables
        - Reservoir level
        - Inflows of previous time periods (according to the lag)
    '''
    #Reservoir level
    reservoir_level = m.addVars(nr, 
                                lb = [r.min for r in valley_chain], 
                                ub = [r.max for r in valley_chain], 
                                obj = 0,
                                vtype=GRB.CONTINUOUS, 
                                name='reservoir_level')
    reservoir_level0 = m.addVars(nr, 
                                lb = 0, 
                                ub = 0, 
                                obj = 0,
                                vtype=GRB.CONTINUOUS, 
                                name='reservoir_level0')
    lag_set = list(range(1,lag+1))
    inflow = m.addVars(nr, lag_set ,lb=-GRB.INFINITY, ub=GRB.INFINITY, obj=0, vtype=GRB.CONTINUOUS, name='inflow')
    inflow0 = m.addVars(nr, lag_set,  lb=0, ub=0, obj=0, vtype=GRB.CONTINUOUS, name='inflow0')
    
    #RHS noise
    innovations = m.addVars(nr,lb=-GRB.INFINITY, ub=GRB.INFINITY,  obj=0, vtype=GRB.CONTINUOUS, name='innovations')
    
    outflow = m.addVars(nr, lb=0, obj=0, vtype=GRB.CONTINUOUS, name='outflow')
    spill = m.addVars(nr, lb=0, obj=0, vtype=GRB.CONTINUOUS, name='spill')
    pour = m.addVars(nr, lb=0, obj=0, vtype=GRB.CONTINUOUS, name='pour')
    generation = m.addVar(lb=0, obj=0, vtype=GRB.CONTINUOUS, name='generation')
    dispatch = m.addVars([(ri,tf)  for (ri,r) in enumerate(valley_chain) for tf in range(0,len(r.turbine.flowknots))],
                        lb=0,ub=1,obj=0,vtype= GRB.CONTINUOUS,name='dispatch')
    if stage == 0:
        for v in reservoir_level0:
            reservoir_level0[v].lb = valley_chain[v].initial
            reservoir_level0[v].ub = valley_chain[v].initial
        for v in inflow0:
            v_lag = v[1] #Lag id
            v_res = v[0] #Reservoir id
            inflow0[v].lb = initial_inflow[-v_lag-1][v_res]
            inflow0[v].ub = initial_inflow[-v_lag-1][v_res]
            inflow[v].lb = initial_inflow[-v_lag][v_res]
            inflow[v].ub = initial_inflow[-v_lag][v_res]
            
    m.update()
            
    in_state = [v.VarName for v in reservoir_level0.values()]
    in_state.extend((v.VarName for v in inflow0.values()))
    out_state = [v.VarName for v in reservoir_level.values()]
    out_state.extend((v.VarName for v in inflow.values()))
    rhs_vars = [v.VarName for v in innovations.values()]
    
    #Constraints
    #AR model for the stage
    R_t = Rmatrix[stage] #For lag 1 only!
    m.addConstrs((inflow[i,1] == sum(R_t[l][i][j]*inflow0[j,l]  for l in lag_set for j in R_t[l][i]) + innovations[i] for i in range(0,len(valley_chain)) ), 'AR_model_%i' %(1))
    m.addConstrs((inflow[i,l] == inflow0[i,l-1]  for l in range(2,lag+1) for i in range(0,len(valley_chain))), 'AR_model')
    
    
    #R_t = Rmatrix[stage] #For lag 1 only!
    #m.addConstrs((inflow[i] - sum(R_t[1][i][j]*inflow0[j]  for j in range(0,len(valley_chain)) if j in R_t[1][i]) == innovations[i]    for i in range(0,len(valley_chain)) ), 'AR1')
    #Balance constraints
    #===========================================================================
    # m.addConstr(reservoir_level[0] ==  reservoir_level0[0] + sum(R_t[l][0][j]*inflow0[j,l]  for l in lag_set for j in R_t[l][0]) + innovations[0] - outflow[0] - spill[0] + pour[0], 'balance[0]')
    # m.addConstrs((reservoir_level[i] ==  reservoir_level0[i] + sum(R_t[l][i][j]*inflow0[j,l]  for l in lag_set for j in R_t[l][i])+ innovations[i] - outflow[i] - spill[i] + pour[i] + outflow[i-1] + spill[i-1] for i in range(1,nr)), 'balance')
    #===========================================================================
         
    m.addConstr(reservoir_level[0] ==  reservoir_level0[0] + inflow[0,1] - outflow[0] - spill[0] + pour[0], 'balance[0]')
    m.addConstrs((reservoir_level[i] ==  reservoir_level0[i] + inflow[i,1] - outflow[i] - spill[i] + pour[i] + outflow[i-1] + spill[i-1] for i in range(1,nr)), 'balance') 
    
    
    #Generation
    m.addConstr(generation==quicksum(r.turbine.powerknots[level] * dispatch[i,level] for (i,r) in enumerate(valley_chain) for level in range(0,len(r.turbine.flowknots))), 'generationCtr')

    # Flow out
    for (i,r) in enumerate(valley_chain):
        m.addConstr(outflow[i] == quicksum(r.turbine.flowknots[level] * dispatch[i, level] for level in range(len(r.turbine.flowknots))), 'outflowCtr[%i]' %(i))
    
    #Dispatched
    for (i,r) in enumerate(valley_chain):
        m.addConstr(quicksum(dispatch[i, level] for level in range(len(r.turbine.flowknots)))<= 1, 'dispatchCtr[%i]' %(i))
    #Objective
    objfun = -prices[stage]*generation + quicksum(0*r.spill_cost*spill[i] for (i,r) in enumerate(valley_chain)) + quicksum(r.spill_cost*pour[i] for (i,r) in enumerate(valley_chain))
    m.setObjective(objfun, GRB.MINIMIZE)
    m.update()
    if stage == -1:
        print(initial_inflow)
        print_model(m)
    return m, in_state, out_state, rhs_vars

def print_model(m):
    for c in m.getConstrs(): print(c.ConstrName, m.getRow(c) , '  ', c.Sense, '  ', c.RHS)
    for v in m.getVars(): print(v.varname, ' '  , v.lb , '  ---  ', v.ub)
    #for v in m.getVars(): print(v)
    
if __name__ == '__main__':
    argv = sys.argv
    positional_args,kwargs = parse_args(argv[1:])
    if 'R' in kwargs:
        nr = kwargs['R']
    if 'T' in kwargs:
        T = kwargs['T']
    if 'max_iter' in kwargs:
        CutSharing.options['max_iter'] = kwargs['max_iter']
        CutSharing.options['lines_freq'] = int(CutSharing.options['max_iter']/10)
    if 'sim_iter' in kwargs:
        CutSharing.options['sim_iter'] = kwargs['sim_iter']
    if 'lag' in kwargs:
        lag = kwargs['lag']
    if 'dro_radius' in kwargs:
        dro_radius = kwargs['dro_radius']
        
    sddp_log.addHandler(logging.FileHandler("HydroAR%i_ESS.log" %(lag), mode='w'))
    hydro_instance = read_instance('hydro_rnd_instance_R10_UD1_T120_LAG1_OUT10K_AR.pkl' , lag = lag)
    
    for nr in [10]:#,10,50,100,500,1000]:
        instance_name = "Hydro_R%i_AR%i_T%i_I%i_ESS" % (nr, lag, T, CutSharing.options['max_iter'])
        Rmatrix = hydro_instance.ar_matrices
        RHSnoise_density = hydro_instance.RHS_noise[0:nr]
        for N_training in [5,10,30]:#[2,3,5,10,20,30]:
            #Reset experiment design stream 
            reset_experiment_desing_gen()
            train_indeces = set(experiment_desing_gen.choice(range(len(RHSnoise_density[0])),size=N_training, replace = False))
            
         
            
            test_indeces = set(range(len(RHSnoise_density[0]))) - train_indeces
            
            assert len(train_indeces.intersection(test_indeces))==0,  'Not disjoint'
            l_train = list(train_indeces)
            l_train.sort()
            RHSnoise = RHSnoise_density[:,l_train]
            dim_p = len(RHSnoise[0]) 
            q_prob = 1/len(RHSnoise[0])
        
       

            initial_inflow = np.array(hydro_instance.inital_inflows)[:,0:nr]
            valley_turbines  = Turbine([50, 60, 70], [55, 65, 70])
            
            
            #For out of sample performance measure
            l_test = list(test_indeces) 
            l_test.sort()
            RHSnoise_oos = RHSnoise_density#[:,l_test]
            valley_chain_oos = [Reservoir(30, 200, 50, valley_turbines, Water_Penalty, x) for x in RHSnoise_oos]
            out_of_sample_rnd_cont = random_builder_out_of_sample(valley_chain_oos)
            
            wasser_valley_chain =[Reservoir(30, 200, 50, valley_turbines, Water_Penalty, x) for x in RHSnoise]
            
            prices = [10+round(5*np.sin(x),2) for x in range(0,T)]
            
            sim_results_com = []
            '''
            Wasserstein DUS Experiment 3
            Uses dat points for origins and continuum support for destinations
            '''
            valley_chain = [Reservoir(30, 200, 50, valley_turbines, Water_Penalty, x) for x in RHSnoise]
            CutSharing.options['multicut'] = True
            instance_name = "Hydro_R%i_AR%i_T%i_I%i_N%iESS" % (nr, lag, T, CutSharing.options['max_iter'], len(valley_chain[0].inflows))
            sim_results = list()
            #for rr in r_lbs:
                
            #for rr in [b*(10**c) for c in [-4,-3,-2,-1,-0,1,2,3] for b in [1,1.5,2,3,4,5,6,7,8,9]]:
            #for rr in [b*(10**c) for c in [1] for b in [1,1.5,2,3,5,9]]:
            rr = dro_radius
            print('Wasserstein Cont r = %10.4e' %(rr))
            #supp_ctrs = [{'innovations[%i]' %(resv):1 for resv in range(nr)} , {'innovations[%i]' %(resv):-1 for resv in range(nr)}]
            #supp_rhs = [RHSnoise_wasswer.sum(axis=0).max(), -(RHSnoise_wasswer.sum(axis=0).min())]             
            supp_ctrs = [{'innovations[%i]' %(resv):1} for resv in range(nr)]
            supp_ctrs.extend(({'innovations[%i]' %(resv):-1}) for resv in range(nr))
            supp_rhs = [RHSnoise[resv].max() for resv in range(nr)] 
            supp_rhs.extend((-RHSnoise[resv].min() for resv in range(nr)))                                                                          
            algo = SDDP(T, model_builder, random_builder, risk_measure = DistRobustWassersteinCont, radius = rr, support_ctrs = supp_ctrs,  support_rhs = supp_rhs)
            algo.run( instance_name=instance_name, dynamic_sampling=True)
            
            sim_result = algo.simulate_policy(CutSharing.options['sim_iter'], out_of_sample_rnd_cont)
            save_path = hydro_path+'/Output/WassersteinCont/%s_WC_%.5f.pickle' %(instance_name,rr)
            write_object_results(save_path, sim_result)
            #sim_results.append(sim_result, sim_result)
            del(algo)
            #plot_sim_results(sim_results, hydro_path+'/Output/WassersteinCont/%s_WassersteinConts.pdf' %(instance_name), N_training)
            #sim_results_com.append(sim_results)
            #plot_metrics_comparison(sim_results_com, hydro_path+'/Output/%s_WassersteinCompX3.pdf' %(instance_name))
            
    

