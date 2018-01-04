'''
Created on Jan 2, 2018

@author: dduque

Generates an instance of the hydro scheduling problem for a chain of reservoirs.

Outputs:
Autoregressive matrix for each time period
'''
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import pickle
from HydroExamples import  hydro_path

def gen_instance(num_reservoirs = 1000, up_stream_dep = 1, T  = 12, lag = 1, num_outcomes = 30 , simulate  = False):
    '''
    Generate a random instance consisting of:
        - Autoregresive matrices (stored as dictionaries)
        - Initial inflow vector
        - Innovations of the autoregressive process
    '''
    np.random.seed(0)
    
    R_matrices = {t:{l:{i:{} for i in range(num_reservoirs)} for l in range(1,lag+1)} for t in range(0,T)}
    for t in range(T):
        for l in range(1,lag+1):
            for i in range(num_reservoirs):
                for j in range(up_stream_dep+1):
                    if i-j>=0:
                        if (t<12):
                            R_matrices[t][l][i][i-j]=np.random.uniform(-0.01/(up_stream_dep+1),2/(up_stream_dep+1))
                        else:
                            R_matrices[t][l][i][i-j]=R_matrices[t-12][l][i][i-j]
    
    inflow_t0 = [np.random.uniform(5,30) for i in range(num_reservoirs)]
    
    RHS_noise = np.zeros(shape=(num_reservoirs,num_outcomes))
    for i in range(num_reservoirs):
        mu_i = mu_ref = np.random.uniform(1,2)
        RHS_noise[i] = np.random.normal(mu_i,1,num_outcomes)
        
    if simulate:
        res_ref  =2
        mu_ref = np.random.uniform(1,2)
        for replica in range(15):
            plot_data = [10]
            inflows = [plot_data[-1]]*num_reservoirs
            for t in range(T):
                innovation  = np.random.normal(mu_ref,1)
                new_inflows = [0]*num_reservoirs
                for i in range(num_reservoirs):
                    for j in range(num_reservoirs):
                        if(j in R_matrices[t][1][i]):
                            new_inflows[i]+=R_matrices[t][1][i][j]*inflows[j]+innovation
                inflows = list(new_inflows)                                         
                plot_data.append(inflows[res_ref])
            plt.plot(plot_data)
        plt.show()
    instance = HydroRndInstance(R_matrices, inflow_t0, RHS_noise)
    return instance

class HydroRndInstance():
    def __init__(self, ar_matrices, initial_inflows, RHS_noise):
        self.ar_matrices = ar_matrices
        self.inital_inflows = initial_inflows
        self.RHS_noise = RHS_noise
        
def read_instance(file_name = 'hydro_rnd_instance_R1000_UD1_T120_LAG1_OUT30.pkl'):
    '''
        Read instance from file and returns a HydroRndInstance object.
    '''
    file_name_path = hydro_path+'/data/'+file_name
    with open(file_name_path, 'rb') as input:
        instance = pickle.load(input)
        return instance
       

if __name__ == '__main__':
    file_name_path = hydro_path+'/data/hydro_rnd_instance_R1000_UD1_T120_LAG1_OUT30.pkl'
    with open(file_name_path, 'wb') as output:
        instance = gen_instance(num_reservoirs=1000, up_stream_dep=1, T=120, lag = 1, num_outcomes=30,  simulate= False)   
        pickle.dump(instance, output, pickle.HIGHEST_PROTOCOL)

