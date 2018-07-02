'''
Created on May 2, 2018

@author: dduque

Module to save and plot simulation results
'''
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
class SimResult():
    '''
    Class the stors the information for a particular instance
    '''
    def __init__(self,instance_params, simulation_upper_bounds):
        self.instance  = instance_params
        self.sims_ub = simulation_upper_bounds
    
    
def plot_lbs2(lbs_list, plot_path):
    assert len(lbs_list) == 2
    dash_styles = [(5, 2),(1, 1),(3, 2),(4, 7)]
    methods_names = ['Empirical', 'Dynamic']  
    x_range = len(max(lbs_list, key= lambda x:len(x)))
    #r = np.arange(x_range)
    f, axarr = plt.subplots(1, 1, figsize=(6, 6), dpi=300)
    min_val = np.inf
    max_val = -np.inf
    for (i,lbs) in enumerate(lbs_list):
        axarr.plot([i for i in range(len(lbs))],lbs, color='black', linestyle='--', dashes=dash_styles[i], label='%s' %(methods_names[i]))
        min_val = np.minimum(min_val, lbs[3])
        max_val = np.maximum(max_val, lbs[-1])
        
    axarr.legend(loc='best', shadow=True, fontsize='small')
    axarr.set_ylim(min_val, max_val)
    #axarr.set_xticks(np.arange(0,x_range+1,1))
    pp = PdfPages(plot_path)
    pp.savefig(f)
    pp.close()
            
def plot_lbs(lb_s, lb_d, N, r_list, plot_path):
    assert len(lb_s) == len(r_list)
    assert len(lb_d) == len(r_list)
    
    dash_styles = [(5, 2),(1, 1),(3, 2),(4, 7)]
    methods_names = ['Empirical', 'Dynamic']  
    
    for (k,r) in enumerate(r_list):
        #r = np.arange(x_range)
        f, axarr = plt.subplots(1, 1, figsize=(6, 6), dpi=300)
        min_val = np.inf
        max_val = -np.inf
        axarr.plot([i for i in range(len(lb_s[k]))],lb_s[k], color='black', linestyle='--', dashes=dash_styles[0], label='%s' %(methods_names[0]))
        axarr.plot([i for i in range(len(lb_d[k]))],lb_d[k], color='black', linestyle='--', dashes=dash_styles[1], label='%s' %(methods_names[1]))
        min_val = np.minimum(lb_s[k][5],lb_s[k][5])
        max_val = np.maximum(lb_s[k][-1],lb_d[k][-1])
            
        axarr.legend(loc='best', shadow=True, fontsize='small')
        axarr.set_ylim(min_val, max_val+0.05*np.abs(max_val-min_val))
        axarr.yaxis.set_minor_locator(MultipleLocator(100))
        axarr.set_xlabel('Iteration')
        axarr.set_ylabel('Lower bound')
        axarr.grid(which='minor', alpha=0.2)
        axarr.grid(which='major', alpha=0.5)
        step = int(len(lb_s[k])/10)
        axarr.set_xticks(np.arange(0,len(lb_s[k]),step))
        complete_file_name = plot_path + 'r_%f.pdf' %(r)
        plt.tight_layout()
        pp = PdfPages(complete_file_name)
        pp.savefig(f)
        pp.close()
        
    '''
    Save to file 
    '''
    experiment_data = pd.DataFrame()
    all_lbs = {'Empirical':lb_s, 'Dynamic':lb_d}
    for (k,r) in enumerate(r_list):
        for m in methods_names:
            col_name = 'lb_%s_r_%f' %(m,r)
            experiment_data[col_name]=all_lbs[m][k]
    
    writer = pd.ExcelWriter('%s.xlsx' %(plot_path))
    sheet_name = 'LBs'
    experiment_data.to_excel(writer,sheet_name)
    writer.save()
    

def plot_sim_results(sim_results, plot_path, N):
    '''
    Plot out-of-sample simulation results
    '''
    r = None
    try:
        r = [sr.instance['risk_measure_params']['radius']  for sr in sim_results]
    except: 
        r = [sr.instance['risk_measure_params']['dro_solver_params']['DUS_radius']  for sr in sim_results]
    mean = [np.mean(sr.sims_ub)  for sr in sim_results]
    median = [np.median(sr.sims_ub)  for sr in sim_results]
    p20 = [np.percentile(sr.sims_ub, q=20)   for sr in sim_results]
    p80 = [np.percentile(sr.sims_ub, q=80)   for sr in sim_results]
    p10 = [np.percentile(sr.sims_ub, q=10)   for sr in sim_results]
    p90 = [np.percentile(sr.sims_ub, q=90)   for sr in sim_results]
    p95 = [np.percentile(sr.sims_ub, q=95)   for sr in sim_results]
    p5 = [np.percentile(sr.sims_ub, q=5)   for sr in sim_results]
    p99 = [np.percentile(sr.sims_ub, q=99)   for sr in sim_results]
    p1 = [np.percentile(sr.sims_ub, q=1)   for sr in sim_results]
    f, axarr = plt.subplots(1, 1, figsize=(6, 6), dpi=300)
    axarr.semilogx(r,mean, color='black', label='Mean')
    axarr.semilogx(r,median, color='black',linestyle='--', label='Median')
    axarr.semilogx(r,p20, color='red', linestyle='--', dashes=(1, 1),label='20-80')
    axarr.semilogx(r,p80, color='red', linestyle='--', dashes=(1, 1))
    axarr.semilogx(r,p10, color='red', linestyle='--', dashes=(3, 1),label='10-90')
    axarr.semilogx(r,p90, color='red', linestyle='--', dashes=(3, 1))
    #===========================================================================
    # axarr.semilogx(r,p5, color='red', linestyle='--', dashes=(7, 3) ,label=' 5 - 95')
    # axarr.semilogx(r,p95, color='red', linestyle='--', dashes=(7, 3))
    # axarr.semilogx(r,p1, color='red', label=' 1 - 99')
    # axarr.semilogx(r,p99, color='red', )
    #===========================================================================
    
    #axarr[1].plot(iterations,test_accuracy, color='b', label='Test acc.')
    #axarr[0].set_ylim([algo_options['opt_tol']*0.95, np.max(train_loss)])
    #axarr[1].set_ylim([0, 1])
    
    axarr.legend(loc='best', shadow=True, fontsize='small')
    #axarr[1].legend(loc='lower right', shadow=True, fontsize='x-large')
    
    # Major ticks every 20, minor ticks every 5
    min_val = -68000#np.round(np.min(p1)-0.01*np.abs(np.min(p1)),-2) - 100
    max_val = -60000#np.round(np.max(p99)+0.01*np.abs(np.max(p1))) + 100
    major_r = 1000#np.abs(max_val-min_val)/10
    minor_r = 100#np.abs(max_val-min_val)/50
    major_ticks = np.arange(min_val, max_val, major_r)
    minor_ticks = np.arange(min_val, max_val, minor_r)
    #===========================================================================
    # axarr.set_ylim(min_val, max_val)
    # axarr.set_yticks(major_ticks)
    # axarr.set_yticks(minor_ticks, minor=True)
    #===========================================================================
    
    # And a corresponding grid
    axarr.grid(which='both')
    
    # Or if you want different settings for the grids:
    axarr.grid(which='minor', alpha=0.2)
    axarr.grid(which='major', alpha=0.5)
    
    axarr.yaxis.set_minor_locator(MultipleLocator(500))
    axarr.set_xlabel('Radius')
    axarr.set_ylabel('Out-of-sample performance')
    axarr.grid(which='minor', alpha=0.2)
    axarr.grid(which='major', alpha=0.5)
    
    plt.tight_layout()
    pp = PdfPages(plot_path)
    pp.savefig(f)
    pp.close()
    
    '''
    Save raw data and percentiles
    '''
    experiment_data = pd.DataFrame()
    experiment_data['Name']=[sr.instance['risk_measure'].__name__ for sr in sim_results]
    #experiment_data['Distance']=[sr.instance['risk_measure_params']['dist_func'].__name__ for sr in sim_results]
    experiment_data['Radius']= r
    experiment_data['Mean']=[np.mean(sr.sims_ub) for sr in sim_results]
    for q in [1,5,10,20,50,80,90,95,99]:
        col_name = 'Percentile%i' %(q)
        experiment_data[col_name]=[np.percentile(sr.sims_ub, q=q, interpolation='lower')   for sr in sim_results]
    n_data_points = len(sim_results[0].sims_ub)
    for i in range(n_data_points):
        col_name = 'Obs%i' %(i)
        experiment_data[col_name]=[sr.sims_ub[i] for sr in sim_results]
    
    writer = pd.ExcelWriter('%s.xlsx' %(plot_path))
    sheet_name = 'Wasserstein%i' %(N)
    try:
        if sim_results[0].instance['risk_measure_params']['dist_func'].__name__ != 'norm_fun':
            sheet_name = 'ChiSquare%i' %(N)
    except:
        pass
    experiment_data.to_excel(writer,sheet_name)
    writer.save()
    
    
def plot_metrics_comparison(sim_results_metrics, plot_path):
    dash_styles = [(5, 1),(1, 1),(3, 2),(4, 7)]
    f, axarr = plt.subplots(1, 1, figsize=(6, 6), dpi=300)
    alg_name = ['Original', 'Extended']
    for (i,sim_results) in enumerate(sim_results_metrics):
        r = None
        try:
            r = [sr.instance['risk_measure_params']['radius']  for sr in sim_results]
        except: 
            r = [sr.instance['risk_measure_params']['dro_solver_params']['DUS_radius']  for sr in sim_results]
        mean = [np.mean(sr.sims_ub)  for sr in sim_results]
        median = [np.median(sr.sims_ub)  for sr in sim_results]
        p20 = [np.percentile(sr.sims_ub, q=20)   for sr in sim_results]
        p80 = [np.percentile(sr.sims_ub, q=80)   for sr in sim_results]
        p10 = [np.percentile(sr.sims_ub, q=10)   for sr in sim_results]
        p90 = [np.percentile(sr.sims_ub, q=90)   for sr in sim_results]
        p95 = [np.percentile(sr.sims_ub, q=95)   for sr in sim_results]
        p5 = [np.percentile(sr.sims_ub, q=5)   for sr in sim_results]
        p99 = [np.percentile(sr.sims_ub, q=99)   for sr in sim_results]
        p1 = [np.percentile(sr.sims_ub, q=1)   for sr in sim_results]
        
        axarr.semilogx(r,mean, color='black', linestyle='--', dashes=dash_styles[i], label='Mean (%s)' %(alg_name[i]))
        #axarr.semilogx(r,median, color='black',linestyle='--', label='Median')
        axarr.semilogx(r,p20, color='red', linestyle='--', dashes=dash_styles[i],label=' 20-80 (%s)' %(alg_name[i]))
        axarr.semilogx(r,p80, color='red', linestyle='--', dashes=dash_styles[i])
        
        #=======================================================================
        # axarr.semilogx(r,p10, color='red', linestyle='--', dashes=(3, 1),label='10-90')
        # axarr.semilogx(r,p90, color='red', linestyle='--', dashes=(3, 1))
        # axarr.semilogx(r,p5, color='red', linestyle='--', dashes=(7, 3) ,label=' 5 - 95')
        # axarr.semilogx(r,p95, color='red', linestyle='--', dashes=(7, 3))
        # axarr.semilogx(r,p1, color='red', label=' 1 - 99')
        # axarr.semilogx(r,p99, color='red', )
        #=======================================================================
        
        #axarr[1].plot(iterations,test_accuracy, color='b', label='Test acc.')
        #axarr[0].set_ylim([algo_options['opt_tol']*0.95, np.max(train_loss)])
        #axarr[1].set_ylim([0, 1])
        
        
        
        
        axarr.legend(loc='best', shadow=True, fontsize='small')
        #axarr[1].legend(loc='lower right', shadow=True, fontsize='x-large')
                
                
                
        # Major ticks every 20, minor ticks every 5
        min_val = -68000#np.round(np.min(p1)-0.01*np.abs(np.min(p1)),-2) - 100
        max_val = -60000#np.round(np.max(p99)+0.01*np.abs(np.max(p1))) + 100
        major_r = 1000#np.abs(max_val-min_val)/10
        minor_r = 100#np.abs(max_val-min_val)/50
        major_ticks = np.arange(min_val, max_val, major_r)
        minor_ticks = np.arange(min_val, max_val, minor_r) 
        
        #=======================================================================
        # axarr.set_yticks(major_ticks)
        # axarr.set_yticks(minor_ticks, minor=True)
        # axarr.set_ylim(min_val, max_val)
        #=======================================================================
        
        # And a corresponding grid
        axarr.grid(which='both')
        
        # Or if you want different settings for the grids:
        axarr.grid(which='minor', alpha=0.2)
        axarr.grid(which='major', alpha=0.5)
        
        axarr.yaxis.set_minor_locator(MultipleLocator(500))
        axarr.set_xlabel('Radius')
        axarr.set_ylabel('Out-of-sample performance')
        axarr.grid(which='minor', alpha=0.2)
        axarr.grid(which='major', alpha=0.5)
    
    plt.tight_layout()
    pp = PdfPages(plot_path)
    pp.savefig(f)
    pp.close()
        