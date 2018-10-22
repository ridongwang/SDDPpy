'''
Created on Jun 13, 2018

@author: dduque
'''
import numpy as np
from CutSharing import logger as sddp_log

def print_model(m):
    '''
    Print out a gurobi model
    '''
    print(m.getConstrs())
    for c in m.getConstrs(): print(c.ConstrName, m.getRow(c) , '  ', c.Sense, '  ', c.RHS)
    for v in m.getVars(): print(v.varname, ' '  , v.lb , '  ---  ', v.ub)
    
    

def report_stats(out_of_sample):
    '''
    Report percentiles and other descriptive statistics of a list of out-of-sample outcomes
    Args:
        out_of_sample (list of float): out-of sample results
    '''
    oos = np.array(out_of_sample)
    p10 = np.percentile(oos,10)
    p90 = np.percentile(oos,90)
    o_median = np.median(oos)
    o_mean = np.mean(oos)
    o_sd = np.std(oos)
    
    sddp_log.info('Mean:   \t%10.2f' %(o_mean))
    sddp_log.info('Median: \t%10.2f' %(o_median))
    sddp_log.info('SD:     \t%10.2f' %(o_sd))
    sddp_log.info('10-90: \t(%10.2f, %10.2f)' %(p10,p90))
