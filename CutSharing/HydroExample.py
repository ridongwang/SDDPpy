'''
Created on Nov 18, 2017

@author: dduque
'''
from gurobipy import *
from CutSharing.SDPP_Alg import SDDP
from CutSharing.RandomnessHandler import *
'''
Instance data
'''
T = 3
nr = 2 #Number of reservoirs
class Turbine():
    def __init__(self, flowknots, powerknots):
        self.flowknots = flowknots
        self.powerknots = powerknots
class Reservoir():
    def __init__(self, minlevel, maxlevel, initial, turbine, s_cost, inflows):
        self.min = minlevel
        self.max = maxlevel
        self.initial = initial
        self.turbine = turbine
        self.spill_cost = s_cost
        self.inflows = inflows 
valley_chain = [
        Reservoir(0, 200, 200, Turbine([50, 60, 70], [55, 65, 70]), 1000, [0, 20, 50]),
        Reservoir(0, 200, 200, Turbine([50, 60, 70], [55, 65, 70]), 1000, [0, 0,  20])]
prices = [1,2,3]


def random_builder():
    rndVectors = []
    for t in range(0,T):
        rv_t = StageRandomVector(t)
        for (i,r) in enumerate(valley_chain):
            rv_t.addRandomElememnt(ele_name='inflow[%i]' %(i), ele_outcomes=r.inflows)    
        
        rndVectors.append(rv_t)
    return rndVectors

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
    outflow = m.addVars(nr, lb=0, obj=0, vtype=GRB.CONTINUOUS, name='outflow')
    inflow = m.addVars(nr, lb=0, ub=0, obj=0, vtype=GRB.CONTINUOUS, name='inflow')
    spill = m.addVars(nr, lb=0, obj=0, vtype=GRB.CONTINUOUS, name='spill')
    generation = m.addVar(lb=0, obj=0, vtype=GRB.CONTINUOUS, name='generation')
    dispatch = m.addVars([(ri,tf)  for (ri,r) in enumerate(valley_chain) for tf in range(0,len(r.turbine.flowknots))],
                        lb=0,ub=1,obj=0,vtype= GRB.CONTINUOUS,name='dispatch')
    if stage == 0:
        for v in reservoir_level0:
            reservoir_level0[v].lb = valley_chain[v].initial
            reservoir_level0[v].ub = valley_chain[v].initial
    m.update()
            
    in_state = [v.VarName for v in reservoir_level0.values()]
    out_state = [v.VarName for v in reservoir_level.values()]
    rhs_vars = [v.VarName for v in inflow.values()]
    
    
    
    #Constraints
    #Balance constraints
    m.addConstr(reservoir_level[0] ==  reservoir_level0[0] + inflow[0] - outflow[0] - spill[0], 'balance[0]')
    m.addConstrs((reservoir_level[i] ==  reservoir_level0[i] + inflow[i] - outflow[i] - spill[i] + outflow[i-1] + spill[i-1]     for i in range(1,nr)), 'balance')
          
    #Generation
    m.addConstr(generation==quicksum(r.turbine.powerknots[level] * dispatch[i,level] for (i,r) in enumerate(valley_chain) for level in range(0,len(r.turbine.flowknots))), 'generationCtr')

    # Flow out
    for (i,r) in enumerate(valley_chain):
        m.addConstr(outflow[i] == quicksum(r.turbine.flowknots[level] * dispatch[i, level] for level in range(len(r.turbine.flowknots))), 'outflowCtr[%i]' %(i))
    
    #Dispatched
    for (i,r) in enumerate(valley_chain):
        m.addConstr(quicksum(dispatch[i, level] for level in range(len(r.turbine.flowknots)))<= 1, 'dispatchCtr[%i]' %(i))
    
    objfun = -prices[stage]*generation + quicksum(r.spill_cost*spill[i] for (i,r) in enumerate(valley_chain))
    m.setObjective(objfun, GRB.MINIMIZE)
    m.update()
    print(m.getObjective())
    print(m.getConstrs())
    return m, in_state, out_state, rhs_vars


    

if __name__ == '__main__':
    algo = SDDP(T, model_builder, random_builder)
    algo.run()
    print(algo)

    
    
    
    
