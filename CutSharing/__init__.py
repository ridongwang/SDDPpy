from gurobipy import *
import numpy as np
import os
import logging


sddp_dir_path = os.path.dirname(os.path.realpath(__file__))
cwd = os.getcwd()


logger = logging.getLogger('SDDP')
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())

'''
Logging settings
'''


'''
Type of passes
'''
FORWARD_PASS = 1
BACKWARD_PASS = 2



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

'''
Cut selector
'''
LAST_CUTS_SELECTOR = 'last_cut_selector'
SLACK_BASED_CUT_SELECTOR = 'slack_based_cut_selector'


'''
Algorithm options
'''
options = {}
options['max_iter'] = 200
options['sim_iter'] = 1000
options['outputlevel']  = 2
options['lines_freq']  = 1
options['n_sample_paths'] = 1
options['grb_threads'] = 1
options['multicut'] = False
options['in_sample_ub'] = 200
options['opt_tol'] = 1E-4
options['dynamic_sampling'] = False
options['max_stage_with_oracle'] = 3
options['max_iters_oracle_ini'] = 10
options['max_cuts_last_cuts_selector'] = 2000
options['slack_cut_selector'] = 1E-4
options['slack_num_iters_cut_selector'] = 300
options['cut_selector'] = None

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
    return options


class not_optimal_sp(Exception):
    def __init__(self, msg):
        super(not_optimal_sp, self).__init__(msg)
        self.msg = msg
    
    
        
        