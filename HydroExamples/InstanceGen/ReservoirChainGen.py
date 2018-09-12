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
                            var = 0.2 if i>num_reservoirs/2 else 0.6
                            R_matrices[t][l][i][i-j]=np.random.normal(0, var)
                            #R_matrices[t][l][i][i-j]=np.random.normal(0.1, (1.0/(lag*up_stream_dep+1)))
                            #R_matrices[t][l][i][i-j]=np.random.uniform(-var,var)
                            #===================================================
                            # if t>0:
                            #     R_matrices[t][l][i][i-j]=np.abs(np.random.normal(0.01, (1.0/(lag*up_stream_dep+1))))
                            #     #R_matrices[t][l][i][i-j]=np.random.uniform(0.0/(up_stream_dep+1)**lag,(1/(up_stream_dep+1))/lag)
                            #     #R_matrices[t][l][i][i-j]=(np.random.normal(0.0, (1.0/(up_stream_dep+1))/lag) + R_matrices[t-1][l][i][i-j])/2.0
                            # else:
                            #     R_matrices[t][l][i][i-j]=np.random.normal(0.01, (1.0/(lag*up_stream_dep+1)))
                            #===================================================
                        else:
                            R_matrices[t][l][i][i-j]=R_matrices[t-12][l][i][i-j]
    
    np.random.seed(1234)
    inflow_t0 = [[np.random.uniform(0,15) for i in range(num_reservoirs)] for l in range(lag+1)] 
    
    print(np.array(inflow_t0)[:,0:5])
    RHS_noise = np.zeros(shape=(num_reservoirs,num_outcomes))
    #mu_s = np.random.uniform(5,10,num_reservoirs)
    mu_s = np.random.uniform(0.5,1.5,num_reservoirs)
    sig_s = np.random.uniform(0.2,1.2,num_reservoirs)
    for i in range(num_reservoirs):
        #RHS_noise[i] = np.sort(np.random.normal(mu_s[i],mu_s[i]/3,num_outcomes))
        #RHS_noise[i] = np.random.normal(mu_s[i],mu_s[i]/3,num_outcomes)
        RHS_noise[i] = np.random.lognormal(mu_s[i],sig_s[i],num_outcomes)
        #RHS_noise[i] = np.sort(np.random.lognormal(mu_s[i],0.5,num_outcomes))
    print(np.max(RHS_noise, 1))
    if simulate:
        plt.figure(1)
        num_reps = 500
        res_ref  = [0,2,4,6,7,9]
        np.random.seed(res_ref)
        mean_res_ref = {rr:np.zeros((T+1)) for rr in res_ref}
        for replica in range(num_reps):
            plot_data = {rr:[inflow_t0[-1][rr]] for rr in res_ref}
            inflows = list(inflow_t0)
            for t in range(T):
                #innovation  = np.random.triangular(-1, mu_ref, 4)
                
                new_inflows = [0]*num_reservoirs
                for l in range(1,lag+1):
                    for i in range(num_reservoirs):
                        for j in range(num_reservoirs):
                            if(j in R_matrices[t][l][i]):
                                new_inflows[i]+=R_matrices[t][l][i][j]*inflows[-l][j]+np.random.choice(RHS_noise[i])
                inflows.append(new_inflows)
                inflows.pop(0)     
                for rr in res_ref:                                   
                    plot_data[rr].append(inflows[-1][rr])
            for (i,rr) in enumerate(res_ref):   
                mean_res_ref[rr] = mean_res_ref[rr] + np.array(plot_data[rr])
                plotpos = int('%i1%i' %(len(res_ref), i+1))
                plt.subplot(plotpos)
                plt.plot(plot_data[rr], alpha = 0.5)
        
        for (i,rr) in enumerate(res_ref):   
            mean_res_ref[rr] = mean_res_ref[rr]/num_reps
            plotpos = int('%i1%i' %(len(res_ref), i+1))
            plt.subplot(plotpos)
            plt.plot(mean_res_ref[rr], linewidth=2, color='black',linestyle='--')
            plt.grid()
        plt.show()
    instance = HydroRndInstance(R_matrices, inflow_t0, RHS_noise)
    return instance

class HydroRndInstance():
    def __init__(self, ar_matrices, initial_inflows, RHS_noise):
        self.ar_matrices = ar_matrices
        self.inital_inflows = initial_inflows
        self.RHS_noise = RHS_noise
        
def read_instance(file_name = 'hydro_rnd_instance_R200_UD1_T120_LAG1_OUT10K_AR.pkl', lag = None):
    '''
        Read instance from file and returns a HydroRndInstance object.
    '''
    file_name_path = hydro_path+'/data/'+file_name
    if lag == None:
        file_name_path = hydro_path+'/data/'+'hydro_rnd_instance_R1000_UD1_T120_LAG%i_OUT30_AR.pkl' %(lag)
    print('hola' , file_name)
    with open(file_name_path, 'rb') as input:
        instance = pickle.load(input)
        return instance
       

if __name__ == '__main__':
    for lag in range(1,2):
        file_name_path = hydro_path+'/data/hydro_rnd_instance_R10_UD1_T120_LAG%i_OUT10K_AR.pkl' %(lag)
        with open(file_name_path, 'wb') as output:
            instance = gen_instance(num_reservoirs=10, up_stream_dep=2, T=24, lag = lag, num_outcomes=10000,  simulate= False)   
            pickle.dump(instance, output, pickle.HIGHEST_PROTOCOL)
    #instance = gen_instance(num_reservoirs=10, up_stream_dep=2, T=24, lag = 1, num_outcomes= 10000,  simulate= True)   
