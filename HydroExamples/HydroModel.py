'''
Created on Sep 12, 2018

@author: dduque

Hydro model functions. Unifies the creation of the model for testing
with different risk measures.
'''
from CutSharing.RandomnessHandler import RandomContainer, StageRandomVector
from CutSharing.RandomManager import reset_experiment_desing_gen, experiment_desing_gen
from CutSharing.SDDP_utils import print_model
from CutSharing import logger as sddp_log, options
from gurobipy import Model, GRB, quicksum
from scipy.spatial import ConvexHull
from Utils.argv_parser import sys,parse_args
from scipy.optimize import linprog
import numpy as np
import logging



'''
Objects for Hydro scheduling examples
'''
import os 
import sys
hydro_path = os.path.dirname(os.path.realpath(__file__))


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
#valley_chain = None
#valley_chain_oos = None
prices = None
Water_Penalty = 1000

def random_builder(valley_chain):
    rc = RandomContainer()
    rndVectors = []
    for t in range(0,T):
        rv_t = StageRandomVector(t)
        rc.append(rv_t)
        for (i,r) in enumerate(valley_chain):
            if t>0:
                re = rv_t.addRandomElememnt('innovations[%i]' %(i), r.inflows[:,t])
            else:
                re = rv_t.addRandomElememnt('innovations[%i]' %(i), [0.0])
            rndVectors.append(rv_t)
    rc.preprocess_randomness()
    return rc




def model_builder(stage, valley_chain):
    '''
    Builds a particular instance of a multistage problem
    '''
    import gurobipy as gb
    m = Model('Hydrovalley_%i' %(stage))
    
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
    thermal = m.addVar(lb=0, obj=0, vtype=GRB.CONTINUOUS, name='thermal_gen')
    thermal_cost = m.addVar(lb=0, obj=0, vtype=GRB.CONTINUOUS, name='thermal_gen')
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
    m.addConstrs((reservoir_level[i] ==  reservoir_level0[i] + inflow[i,1] - outflow[i] - spill[i] + pour[i] + outflow[i-1] + spill[i-1] for i in range(1,nr) if i not in [0,5]), 'balance') 
    m.addConstr(reservoir_level[5] ==  reservoir_level0[5] + inflow[5,1] - outflow[5] - spill[5] + pour[5], 'balance[5]')
    
    
    #Hydro generation
    m.addConstr(generation==quicksum(r.turbine.powerknots[level] * dispatch[i,level] for (i,r) in enumerate(valley_chain) for level in range(0,len(r.turbine.flowknots))), 'generationCtr')
    #Demand
    m.addConstr(generation+thermal>=400)
    # Flow out
    for (i,r) in enumerate(valley_chain):
        m.addConstr(outflow[i] == quicksum(r.turbine.flowknots[level] * dispatch[i, level] for level in range(len(r.turbine.flowknots))), 'outflowCtr[%i]' %(i))
    
    #Dispatched
    for (i,r) in enumerate(valley_chain):
        m.addConstr(quicksum(dispatch[i, level] for level in range(len(r.turbine.flowknots)))<= 1, 'dispatchCtr[%i]' %(i))
    #Thermal cost ctr
    m.addConstr(thermal_cost >= 10*thermal)
    m.addConstr(thermal_cost >= 20*thermal-500)
    m.addConstr(thermal_cost >= 30*thermal-2500)
    
    #Objective 
    objfun = -prices[stage]*generation + quicksum(0*r.spill_cost*spill[i] for (i,r) in enumerate(valley_chain)) + quicksum(r.spill_cost*pour[i] for (i,r) in enumerate(valley_chain)) + thermal_cost
    m.setObjective(objfun, GRB.MINIMIZE)
    m.update()
    if stage == -1:
        print(initial_inflow)
        print_model(m)
    return m, in_state, out_state, rhs_vars


def generate_extra_data(rvs, n, method='cvx_hull', realizations=None):
    '''
    Generates additional samples of the random vectors.
    Random vectors are assumed to be in R^d
    Args:
        rvs (ndarray): array of random vectors (Nxd)
        n (int): number of samples to generate
    Returns
        gen_data (ndarray): new data points (nxd)
    '''
    
    if method == 'cvx_hull':
        #Get convex hull first to generate inside if dimensions allow
        ch_points = rvs
        if len(rvs)>len(rvs[0]):
            ch=ConvexHull(rvs)
            ch_points = rvs[ch.vertices]
        #Generate exponential numbers and then normalize to perform the cvx combinations
        cvx_com = experiment_desing_gen.exponential(1,size=(n,len(ch_points)))
        cvx_com = (cvx_com.transpose()/np.sum(cvx_com,1)).transpose()
        #Generated points
        gen = cvx_com.dot(ch_points)
        return gen
    elif method == 'cvx_hull2':
        assert len(rvs)>len(rvs[0]), "Need at least %i points" %(len(rvs[0])+1)
        points_out =[]
        ch=ConvexHull(rvs)
        ch_points = rvs[ch.vertices]
        n_ch = len(ch_points)
        A = np.vstack((ch_points.transpose(), np.ones(n_ch)))
        for b in realizations:
           
            b_e =np.reshape(np.hstack((b,1)),(len(b)+1,1))
            sol = linprog(c=np.zeros(n_ch), A_eq=A, b_eq = b_e, bounds=[(0,1) for _ in range(n_ch)])
            if sol.status == 2 : #infeasible:
                points_out.append(b)
            else:
                assert sol.status == 0, 'Linear system is not infeasible nor optimal.'
                print("interior point")
                pass #Point is inside
            
            if len(points_out)==n:
                break
        return np.array(points_out)
    elif method == 'box':
        
        #cube operations
        d = len(rvs[0])
        min_vals = rvs.min(0)
        max_vals = rvs.max(0)
        #u = experiment_desing_gen.uniform(size=(n,d))
        gen = np.zeros((n,d))
        
        #convex hull operations
        ch=ConvexHull(rvs)
        ch_points = rvs[ch.vertices]
        n_ch = len(ch_points)
        A = np.vstack((ch_points.transpose(), np.ones(n_ch)))
        
        #generate outside the cupe but inside the box
        i = 0
        while i < n:
            u = experiment_desing_gen.uniform(size=(d))
            new_p = u*min_vals + (1-u)*max_vals
            #check if it is insde the cvx hull
            b_e =np.reshape(np.hstack((new_p,1)),(len(new_p)+1,1))
            sol = linprog(c=np.zeros(n_ch), A_eq=A, b_eq = b_e, bounds=[(0,1) for _ in range(n_ch)])
            if sol.status == 2 : #infeasible:
                gen[i] = np.copy(new_p)
                i+=1
            else:
                assert sol.status == 0, 'Linear system is not infeasible nor optimal.'
                print("interior point")
                pass #Point is inside
            
        #=======================================================================
        # for i in range(n):
        #     for k in range(d):
        #         gen[i,k] = u[i,k]*min_vals[k] + (1-u[i,k])*max_vals[k]
        #=======================================================================
        
        return gen

def load_hydro_data(approach, dus_type):
    global T 
    global nr 
    global lag 
    global dro_radius 
    global Rmatrix 
    global RHSnoise 
    global initial_inflow
    global prices 
    argv = sys.argv
    DW_extended = 1
    DW_sampling = None
    positional_args,kwargs = parse_args(argv[1:])
    if 'R' in kwargs:
        nr = kwargs['R']
    if 'T' in kwargs:
        T = kwargs['T']
    if 'lag' in kwargs:
        lag = kwargs['lag']
    if 'dro_r' in kwargs:
        dro_radius = kwargs['dro_r']
    if 'N' in kwargs:
        N = kwargs['N']
    if 'DW_extended' in kwargs:
        #N*DW_extended would be the number of oracles
        DW_extended = kwargs['DW_extended']
    if 'DW_sampling' in kwargs:
        DW_sampling = kwargs['DW_sampling']
        
        
    from InstanceGen.ReservoirChainGen import read_instance
    prices = [10+round(5*np.sin(0.5*(x-2)),2) for x in range(0,T)]
    print(prices)
    hydro_instance = read_instance('hydro_rnd_instance_R30_UD1_T24_LAG1_OUT10K_AR1.pkl' , lag = lag)
    Rmatrix = hydro_instance.ar_matrices
    RHSnoise_density = hydro_instance.RHS_noise[0:nr, : , 0:T] #Total of 10,000 samples
    initial_inflow = np.array(hydro_instance.inital_inflows)[:,0:nr]
    valley_turbines  = Turbine([50, 60, 70], [55, 65, 70])
    
    N_data = N
    #Reset experiment design stream 
    reset_experiment_desing_gen()
    
    #For out of sample performance measure
    test_indeces = set(experiment_desing_gen.choice(range(len(RHSnoise_density[0])),size=9000, replace = False))
    l_test = list(test_indeces) 
    l_test.sort()
    RHSnoise_oos = RHSnoise_density[:,l_test]
    valley_chain_oos = [Reservoir(30, 200, 50, valley_turbines, Water_Penalty, x) for x in RHSnoise_oos]
    
    #Train indices for Wasserstein distance
    available_indices = set(range(len(RHSnoise_density[0]))) - test_indeces
    available_indices = np.array(list(available_indices))
    #data_indeces=set(experiment_desing_gen.choice(list(available_indices), size=N_data, replace=False))
    data_indeces = available_indices[0:N_data]
    RHSnoise_data = RHSnoise_density[:,data_indeces]
    valley_chain_data = [Reservoir(30, 200, 50, valley_turbines, Water_Penalty, x) for x in RHSnoise_data]
    print(data_indeces)
    if DW_extended > 1 and dus_type=='DW':
        #Generate additional data points from the data
        if DW_sampling==None or DW_sampling=='none'  or DW_sampling=='None':
            #available_indices = set(available_indices) - set(data_indeces)
            N_wasserstein = N_data * DW_extended 
            #train_indeces = set(experiment_desing_gen.choice(list(available_indices), size=N_wasserstein, replace=False)) 
            train_indeces = available_indices[0:N_wasserstein]
            assert set(data_indeces).issubset(train_indeces)
            N_training = len(train_indeces)
            l_train = list(train_indeces)
            l_train.sort()
            RHSnoise = RHSnoise_density[:,l_train]
        else:
            N_wasserstein = N_data * DW_extended - N_data
            available_indices = available_indices[N_data:]
            avail_realizations = RHSnoise_density[:,available_indices]
            gen_data = generate_extra_data(RHSnoise_data.transpose(), N_wasserstein, method=DW_sampling, realizations=avail_realizations.transpose())
            RHSnoise = np.hstack((RHSnoise_data,gen_data.transpose()))
            N_training = len(RHSnoise[0])
    else:
        train_indeces = data_indeces
        N_training = N_data
        RHSnoise = np.copy(RHSnoise_data)    

    rr = dro_radius
    cut_type = 'MC' if options['multicut'] else 'SC'
    sampling_type = 'DS' if options['dynamic_sampling']  else 'ES'
    instance_name = None
    if approach == "SP":
        instance_name = "Hydro_R%i_AR%i_T%i_N%i_%i_I%i_Time%i_%s_%s_%s" % (nr, lag, T, N_data, N_training,  options['max_iter'],options['max_time'], approach, cut_type,  sampling_type)
    else:
        instance_name = "Hydro_R%i_AR%i_T%i_N%i_%i_I%i_Time%i_%s_%s_%s_%s_%.7f_%s" % (nr, lag, T, N_data, N_training,  options['max_iter'],options['max_time'], dus_type, approach, cut_type,  sampling_type , rr, DW_sampling)
    sddp_log.addHandler(logging.FileHandler(hydro_path+"/Output/log/%s.log" %(instance_name), mode='w'))
    print(RHSnoise)
    valley_chain = [Reservoir(30, 200, 50, valley_turbines, Water_Penalty, x) for x in RHSnoise]

    def rnd_builder_n_train():
        return random_builder(valley_chain)
    
    def model_builder_n_tran(stage):
        return model_builder(stage, valley_chain)
    
    rnd_container_oos = random_builder(valley_chain_oos)
    rnd_container_data = random_builder(valley_chain_data)
    
    
    return T, model_builder_n_tran, rnd_builder_n_train, rnd_container_data, rnd_container_oos, rr, instance_name
    
