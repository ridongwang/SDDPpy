from gurobipy import *
import numpy as np




'''
Status codes 
'''
SP_LOADED = 'Loaded'
SP_OPTIMAL = 'sp_optimal'
SP_INFEASIBLE = 'sp_infeasible'
SP_UNKNOWN = 'sp_unknown'

'''
tolerances
'''
ZERO_TOL = 1E-8
SDDP_OPT_TOL = 1E-3




def gurobiStatusCodeToStr( intstatus ):
    if intstatus == 1:
        return SP_LOADED
    elif intstatus ==2:
        return SP_OPTIMAL
    elif intstatus in [3,4]:
        return SP_INFEASIBLE
    else:
        return SP_UNKNOWN


def alg_options():
    options = {}
    options['max_iter'] = 100
    options['outputlevel']  = 2
    options['n_sample_paths'] = 20
    return options


class not_optimal_sp(Exception):
    def __init__(self, msg):
        super(not_optimal_sp, self).__init__(msg)
        self.msg = msg
    
    
        
        