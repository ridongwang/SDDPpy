from gurobipy import *
import numpy as np




'''
Status codes 
'''
SP_LOADED = 'Loaded'
SP_OPTIMAL = 'sp_optimal'
SP_INFEASIBLE = 'sp_infeasible'
SP_UNKNOWN = 'sp_unknown'


ZERO_TOL = 1E-8
SDDP_OPT_TOL = 1E-3




def gurobiStatusCodeToStr( intstatus ):
    if intstatus == 1:
        return SP_LOADED
    elif intstatus ==2:
        return SP_OPTIMAL
    elif intstatus == 3:
        return SP_INFEASIBLE
    elif intstatus == 4:
        return SP_UNKNOWN


def alg_options():
    options = {}
    options['max_iter'] = 1500
    options['outputlevel']  = 2
    return options


class not_optimal_sp(Exception):
    def __init__(self, msg):
        super(not_optimal_sp, self).__init__(msg)
        self.msg = msg
    
    
        
        