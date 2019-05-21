'''
Created on Jan 2, 2018

@author: dduque

Generates an instance of the hydro scheduling problem for a chain of reservoirs.

Outputs:
Autoregressive matrix for each time period
'''
if __name__ == '__main__':
    from os import path
    import sys
    sys.path.append(path.abspath('/Users/dduque/Dropbox/WORKSPACE/SDDP'))
    sys.path.append(path.abspath('/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples'))
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import pickle
from HydroModel import hydro_path

def gen_instance(num_reservoirs = 1000, up_stream_dep = 1, T  = 12, lag = 1, num_outcomes = 30 , simulate  = False):
    '''
    Generate a random instance consisting of:
        - Autoregresive matrices (stored as dictionaries)
        - Initial inflow vector (matrix for lag > 0)
        - Innovations of the autoregressive process
    '''
    np.random.seed(0)
    season = 12
    R_matrices = {t:{l:{i:{} for i in range(num_reservoirs)} for l in range(1,lag+1)} for t in range(0,T)}
    for t in range(T):
        for l in range(1,lag+1):
            for i in range(num_reservoirs):
                for j in range(up_stream_dep+1):
                    if i-j>=0:
                        if (t<season):
                            #var = 0.2 if i>num_reservoirs/2 else 0.6
                            #R_matrices[t][l][i][i-j]=np.random.normal(0, var) #for nr=10  experiments
                            R_matrices[t][l][i][i-j]=np.random.normal(0.4, 0.4) #for nr=30  experiments
                            #R_matrices[t][l][i][i-j]=np.random.normal(0.1, (1.0/(lag*up_stream_dep+1)))
                            #R_matrices[t][l][i][i-j]=np.random.uniform(-var,var)
                            #R_matrices[t][l][i][i-j]=np.random.uniform(-1/(up_stream_dep+lag),1/(up_stream_dep+lag)) #for nr=100
                            #===================================================
                            # if t>0:
                            #     R_matrices[t][l][i][i-j]=np.abs(np.random.normal(0.01, (1.0/(lag*up_stream_dep+1))))
                            #     #R_matrices[t][l][i][i-j]=np.random.uniform(0.0/(up_stream_dep+1)**lag,(1/(up_stream_dep+1))/lag)
                            #     #R_matrices[t][l][i][i-j]=(np.random.normal(0.0, (1.0/(up_stream_dep+1))/lag) + R_matrices[t-1][l][i][i-j])/2.0
                            # else:
                            #     R_matrices[t][l][i][i-j]=np.random.normal(0.01, (1.0/(lag*up_stream_dep+1)))
                            #===================================================
                        else:
                            R_matrices[t][l][i][i-j]=R_matrices[t-season][l][i][i-j]
    print(R_matrices[2][1])
    np.random.seed(1234)
    inflow_t0 = [[np.random.uniform(5,40) for i in range(num_reservoirs)] for l in range(lag+1)] 
    
    print(np.array(inflow_t0)[:,0:5])
    #===========================================================================
    # RHS_noise = np.zeros(shape=(num_reservoirs,num_outcomes))
    # #mu_s = np.random.uniform(5,10,num_reservoirs)
    # mu_s = np.random.uniform(0.5,1.5,num_reservoirs)
    # sig_s = np.random.uniform(0.2,1.2,num_reservoirs)
    # for i in range(num_reservoirs):
    #     #RHS_noise[i] = np.sort(np.random.normal(mu_s[i],mu_s[i]/3,num_outcomes))
    #     #RHS_noise[i] = np.random.uniform(-5,10, num_outcomes)
    #     #RHS_noise[i] = np.random.normal(8,np.log(num_reservoirs-i)+1, num_outcomes)
    #     RHS_noise[i] = np.random.lognormal(mu_s[i],sig_s[i],num_outcomes) #nr 10 and nr 100
    #     #RHS_noise[i] = np.sort(np.random.lognormal(mu_s[i],0.5,num_outcomes))
    # print(np.max(RHS_noise, 1))
    #===========================================================================
    RHS_noise = np.zeros(shape=(num_reservoirs,num_outcomes,T))
    for t in range(T):
        mean_t =  np.minimum(np.array([1.5 - round(0.1 * np.sin(0.5 * (t - 2)), 2) for i in range(num_reservoirs)]),1.5)
        sig_t = np.array([1 + round(0.3 * np.sin(0.5 * (t - 2)), 2) for i in range(num_reservoirs)])
        
        
        print(mean_t[0], '  ' , sig_t[0])
        mu_s = np.random.uniform(mean_t , mean_t, num_reservoirs)
        sig_s = np.random.uniform(sig_t*0.5 , sig_t, num_reservoirs)
        #loc_s = np.exp(mu_s+0.5*sig_s**2)
        for i in range(num_reservoirs):
            if t<season:
                #RHS_noise[i,:,t] = np.random.lognormal(mu_s[i],sig_s[i],num_outcomes) #nr 10 and nr 100
                mu_i = np.exp(mu_s[i]+0.5*sig_s[i]**2)
                var_i = mu_i*sig_s[i]
                RHS_noise[i,:,t] = np.abs(np.random.normal(mu_i, var_i,num_outcomes)) #nr 10 and nr 100
            else:
                RHS_noise[i,:,t] = RHS_noise[i,:,t-season]
        
    
    if simulate:
        simulate_AR_model(R_matrices, inflow_t0, RHS_noise, T, num_reservoirs, lag)
        
    instance = HydroRndInstance(R_matrices, inflow_t0, RHS_noise)
    return instance



def simulate_AR_model(R_matrices, inflow_t0, RHS_noise, T, nr, lag):
    prices = [10+round(5*np.sin(0.5*(x-2)),2) for x in range(0,T)]
    num_reservoirs = nr
    plt.figure(1)
    num_reps = 200
    res_ref  = [0, 5, 8]
    np.random.seed(res_ref)
    mean_res_ref = {rr:np.zeros((T)) for rr in res_ref}
    for replica in range(num_reps):
        
        plot_data = {rr:[inflow_t0[-1][rr]] for rr in res_ref}
        inflows = list(inflow_t0)
        for t in range(1,T):
            #innovation  = np.random.triangular(-1, mu_ref, 4)
            outcome = np.random.randint(len(RHS_noise[0,:,t]))
            new_inflows = [0]*num_reservoirs
            for l in range(1,lag+1):
                for i in range(num_reservoirs):
                    new_inflows[i]+= RHS_noise[i,outcome,t]
                    for j in range(num_reservoirs):
                        if(j in R_matrices[t][l][i]):
                            new_inflows[i]+=R_matrices[t][l][i][j]*inflows[-l][j]
            inflows.append(new_inflows)
            inflows.pop(0)     
            for rr in res_ref:                                   
                plot_data[rr].append(inflows[-1][rr])
        for (i,rr) in enumerate(res_ref):   
            mean_res_ref[rr] = mean_res_ref[rr] + np.array(plot_data[rr])
            plotpos = int('%i1%i' %(len(res_ref), i+1))
            plt.subplot(plotpos)
            plt.plot(plot_data[rr], alpha = 0.5)
        data_replica = np.array([plot_data[r1] for r1 in res_ref])
    for (i,rr) in enumerate(res_ref):   
        mean_res_ref[rr] = mean_res_ref[rr]/num_reps
        plotpos = int('%i1%i' %(len(res_ref), i+1))
        plt.subplot(plotpos)
        plt.plot(mean_res_ref[rr], linewidth=2, color='black',linestyle='--')
        plt.plot(prices, linewidth=2, color='red',linestyle='--')
        plt.grid()
    plt.show()
        
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
    with open(file_name_path, 'rb') as input:
        instance = pickle.load(input)
        return instance
       

if __name__ == '__main__':
    #===========================================================================
    # hydro_instance = read_instance('hydro_rnd_instance_R30_UD1_T24_LAG1_OUT10K_AR1.pkl' , lag=1)
    # matrix = hydro_instance.ar_matrices
    # RHSnoise_density = hydro_instance.RHS_noise
    # inflow_t0 = hydro_instance.inital_inflows
    # simulate_AR_model(matrix, inflow_t0, RHSnoise_density, 24, 10, 1)
    #===========================================================================
    
    nr = 30
    ud = 1
    for lag in []:#range(1,2):
        file_name_path = hydro_path+'/data/hydro_rnd_instance_R%i_UD%i_T24_LAG%i_OUT10K_AR1.pkl' %(nr,ud,lag)
        print(file_name_path)
        with open(file_name_path, 'wb') as output:
            instance = gen_instance(num_reservoirs=nr, up_stream_dep=ud, T=24, lag = lag, num_outcomes=10000,  simulate=False)   
            pickle.dump(instance, output, pickle.HIGHEST_PROTOCOL)
    instance = gen_instance(num_reservoirs=nr, up_stream_dep=ud, T=24, lag = 1, num_outcomes= 10000,  simulate= True)   
