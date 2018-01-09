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
        - Initial inflow vector (matrix for lag > 0)
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
                            #R_matrices[t][l][i][i-j]=np.random.uniform(0.0/(up_stream_dep+1)**lag,(1/(up_stream_dep+1))/lag)
                            if t>0:
                                R_matrices[t][l][i][i-j]=(np.random.normal(0.1, (1.0/(up_stream_dep+1))/lag) + R_matrices[t-1][l][i][i-j])/2.0
                            else:
                                R_matrices[t][l][i][i-j]=np.random.normal(0.1, (1.0/(up_stream_dep+1))/lag)
                        else:
                            R_matrices[t][l][i][i-j]=R_matrices[t-12][l][i][i-j]
    
    inflow_t0 = [[np.random.uniform(0,30) for i in range(num_reservoirs)] for l in range(lag)] 
    
    RHS_noise = np.zeros(shape=(num_reservoirs,num_outcomes))
    mu_s = np.random.uniform(1,1.5,num_reservoirs)
    for i in range(num_reservoirs):
        RHS_noise[i] = np.sort(np.random.normal(mu_s[i],1,num_outcomes))
        
    if simulate:
        num_reps = 50
        res_ref  = 86
        np.random.seed(res_ref)
        mean_res_ref = np.zeros((T+1))
        for replica in range(num_reps):
            plot_data = [inflow_t0[-1][res_ref]]
            inflows = list(inflow_t0)
            for t in range(T):
                innovation  = np.random.choice(RHS_noise[i])
                #innovation  = np.random.triangular(-1, mu_ref, 4)
                
                new_inflows = [0]*num_reservoirs
                for l in range(1,lag+1):
                    for i in range(num_reservoirs):
                        for j in range(num_reservoirs):
                            if(j in R_matrices[t][1][i]):
                                new_inflows[i]+=R_matrices[t][l][i][j]*inflows[-l][j]+innovation
                inflows.append(new_inflows)
                inflows.pop(0)                                       
                plot_data.append(inflows[-1][res_ref])
            mean_res_ref = mean_res_ref + np.array(plot_data)
            plt.plot(plot_data, alpha = 0.5)
        mean_res_ref = mean_res_ref/num_reps
        plt.plot(mean_res_ref, linewidth=2, color='black',linestyle='--')
        
        plt.show()
    instance = HydroRndInstance(R_matrices, inflow_t0, RHS_noise)
    return instance

class HydroRndInstance():
    def __init__(self, ar_matrices, initial_inflows, RHS_noise):
        self.ar_matrices = ar_matrices
        self.inital_inflows = initial_inflows
        self.RHS_noise = RHS_noise
        
def read_instance(file_name = 'hydro_rnd_instance_R1000_UD1_T120_LAG1_OUT30_V2.pkl', lag = None):
    '''
        Read instance from file and returns a HydroRndInstance object.
    '''
    file_name_path = hydro_path+'/data/'+file_name
    if lag != None:
        file_name_path = hydro_path+'/data/'+'hydro_rnd_instance_R1000_UD1_T120_LAG%i_OUT30_AR.pkl' %(lag)
    with open(file_name_path, 'rb') as input:
        instance = pickle.load(input)
        return instance
       

if __name__ == '__main__':
    for lag in range(1,7):
        file_name_path = hydro_path+'/data/hydro_rnd_instance_R1000_UD1_T120_LAG%i_OUT30_AR.pkl' %(lag)
        with open(file_name_path, 'wb') as output:
            instance = gen_instance(num_reservoirs=1000, up_stream_dep=1, T=120, lag = lag, num_outcomes=30,  simulate= False)   
            pickle.dump(instance, output, pickle.HIGHEST_PROTOCOL)
    #instance = gen_instance(num_reservoirs=100, up_stream_dep=1, T=24, lag =3, num_outcomes=30,  simulate= True)   
