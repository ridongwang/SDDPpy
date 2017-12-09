'''
Created on Nov 18, 2017

@author: dduque
'''
import CutSharing

'''
Instance data
'''
import csv
from gurobipy import *
import numpy as np
from CutSharing.RandomnessHandler import RandomContainer,StageRandomVector,AR1_depedency
from CutSharing.SDPP_Alg import SDDP
from HydroExamples import Reservoir, Turbine 
from ScenarioReduction.ScenarioObjects import Scenario
from ScenarioReduction.SRAlg import kmeans_reduce_set
T = 10
m = 50
AR1Matrix = [] 
Rmatrix = []
with open('../TimeSeries/AR1Matrix%i.csv' %(m), 'r') as f:
    reader = csv.reader(f)
    AR1Matrix = list(reader)
    AR1Matrix.pop(0)
    for (i,x) in enumerate(AR1Matrix):
        x.pop(0)
        Rmatrix.append({})  
        for (j,val) in enumerate(x):
            x[j] = float(val)  
            Rmatrix[-1]['inflow[%i]' %(j)] = x[j]
RHSnoise = [] 
with open('../TimeSeries/RHSnoise%i.csv' %(m), 'r') as f:
    reader = csv.reader(f)
    RHSnoise = list(reader)
    RHSnoise.pop(0)
    for (i,x) in enumerate(RHSnoise):
        x.pop(0) 
        for (j,val) in enumerate(x):
            x[j] = float(val)
            
valley_chain = [
        Reservoir(0, 500, 200, Turbine([50, 60, 70], [55, 65, 70]), 1000, [x[0],x[15],x[-1]]) for x in RHSnoise
        ]
valley_chain = valley_chain[0:2]
nr = len(valley_chain) #Number of reservoirs


prices = [np.minimum(x,5) for x in range(0,T)]

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
    return rc

def model_builder(stage):
    '''
    Builds a particular instance of a multistage problem
    '''
    m = Model('Hydrovalley')

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
    inflow = m.addVars(nr,lb=-GRB.INFINITY, ub=GRB.INFINITY, obj=0, vtype=GRB.CONTINUOUS, name='inflow')
    inflow0 = m.addVars(nr, lb=0, ub=0, obj=0, vtype=GRB.CONTINUOUS, name='inflow0')
    
    
    outflow = m.addVars(nr, lb=0, obj=0, vtype=GRB.CONTINUOUS, name='outflow')
    spill = m.addVars(nr, lb=0, obj=0, vtype=GRB.CONTINUOUS, name='spill')
    innovations = m.addVars(nr,  obj=0, vtype=GRB.CONTINUOUS, name='innovations')
    pour = m.addVars(nr, lb=0, obj=0, vtype=GRB.CONTINUOUS, name='pour')
    generation = m.addVar(lb=0, obj=0, vtype=GRB.CONTINUOUS, name='generation')
    dispatch = m.addVars([(ri,tf)  for (ri,r) in enumerate(valley_chain) for tf in range(0,len(r.turbine.flowknots))],
                        lb=0,ub=1,obj=0,vtype= GRB.CONTINUOUS,name='dispatch')
    if stage == 0:
        for v in reservoir_level0:
            reservoir_level0[v].lb = valley_chain[v].initial
            reservoir_level0[v].ub = valley_chain[v].initial
        for v in inflow0:
            inflow0[v].lb = 0
            inflow0[v].ub = 0
            inflow[v].lb = 0
            inflow[v].ub = 0
            
    m.update()
            
    in_state = [v.VarName for v in reservoir_level0.values()]
    in_state.extend((v.VarName for v in inflow0.values()))
    out_state = [v.VarName for v in reservoir_level.values()]
    out_state.extend((v.VarName for v in inflow.values()))
    rhs_vars = [v.VarName for v in innovations.values()]
    
    #Constraints
    #AR1 model as a constraint
    m.addConstrs((inflow[i] - sum(AR1Matrix[i][j]*inflow0[j]  for j in range(0,len(valley_chain))) == innovations[i]    for i in range(0,len(valley_chain)) ), 'AR1')
    #Balance constraints
    m.addConstr(reservoir_level[0] ==  reservoir_level0[0] + inflow[0] - outflow[0] - spill[0] + pour[0], 'balance[0]')
    m.addConstrs((reservoir_level[i] ==  reservoir_level0[i] + inflow[i] - outflow[i] - spill[i] + pour[i] + outflow[i-1] + spill[i-1]     for i in range(1,nr)), 'balance')
          
    #Generation
    m.addConstr(generation==quicksum(r.turbine.powerknots[level] * dispatch[i,level] for (i,r) in enumerate(valley_chain) for level in range(0,len(r.turbine.flowknots))), 'generationCtr')

    # Flow out
    for (i,r) in enumerate(valley_chain):
        m.addConstr(outflow[i] == quicksum(r.turbine.flowknots[level] * dispatch[i, level] for level in range(len(r.turbine.flowknots))), 'outflowCtr[%i]' %(i))
    
    #Dispatched
    for (i,r) in enumerate(valley_chain):
        m.addConstr(quicksum(dispatch[i, level] for level in range(len(r.turbine.flowknots)))<= 1, 'dispatchCtr[%i]' %(i))
    #Objective
    objfun = -prices[stage]*generation + quicksum(r.spill_cost*spill[i] for (i,r) in enumerate(valley_chain)) + quicksum(r.spill_cost*pour[i] for (i,r) in enumerate(valley_chain))
    m.setObjective(objfun, GRB.MINIMIZE)
    m.update()

    return m, in_state, out_state, rhs_vars

if __name__ == '__main__':
    opt = CutSharing.alg_options() 
    print('Scenario reduction approach')
    algo = SDDP(T, model_builder, random_builder)
    print('Begin enumeration')
    all_sce = algo.random_container.enumerate_scenarios()
    names = [n for n in all_sce[0][0]]
    names.sort()
    dim = T*len(names)
    scenario_list = []
    N = len(all_sce)
    for (i,sp) in enumerate(all_sce):
        values = np.array([x[nv] for x in sp for nv in names])
        scenario_list.append(Scenario(i, values, 1.0/N))
    print('Begin reduction: from %i to %i' %(N, opt['max_iter']))
    reduced_set = kmeans_reduce_set(scenario_list, opt['max_iter'] )
    reduced_set_sam_paths = [all_sce[s.id] for s in reduced_set.values()]
    print('Begin SDDP')
    lbs_sce_red = algo.run(pre_sample_paths = reduced_set_sam_paths)
    print('\n\n================================================')
    
    print('Sampling approach')
    print('Begin SDDP')
    algo = SDDP(T, model_builder, random_builder)
    lbs_samlping  = algo.run()
    
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    f, axarr = plt.subplots(1,1,figsize=(10, 7), dpi=300)
    axarr.plot(range(1,opt['max_iter']+1),lbs_sce_red, color='b', label='LB Scenario reduction')
    axarr.plot(range(1,opt['max_iter']+1),lbs_samlping, color='r', label='LB Monte Carlo Sampling')
    axarr.set_ylim([-5000, -4.00e+03])
    #-4.02875e+03 -4.03031e+03
    legend = axarr.legend(loc='lower right', shadow=True, fontsize='x-large')
    pp = PdfPages('../ScenarioReduction/lbs_%i.pdf' %(opt['max_iter']))
    pp.savefig(f)
    pp.close()
    #plt.show()
    
    
