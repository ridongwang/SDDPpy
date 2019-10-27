'''
Created on Sep 13, 2018

@author: dduque
'''
#===============================================================================
# import os
# import sys
# #Include all modules of the project when running on shell
# this_path = os.path.dirname(os.path.abspath(__file__))
# parent_path= os.path.abspath(os.path.join(this_path, os.pardir))
# cwd = os.getcwd()
# print(__file__, this_path, parent_path, cwd, os.path.abspath(os.pardir))
# sys.path.append(parent_path)
#===============================================================================
from __init__ import import_SDDP, dro_radii
#from HydroExamples import dro_radii, import_SDDP
import_SDDP()

from CutSharing import options, LAST_CUTS_SELECTOR, load_algorithm_options,\
    SLACK_BASED_CUT_SELECTOR
from CutSharing.SDDP_Alg import SDDP
from CutSharing.RiskMeasures import DistRobustWasserstein
from Utils.file_savers import write_object_results
from HydroModel import load_hydro_data, hydro_path
from InstanceGen.ReservoirChainGen import read_instance, HydroRndInstance  # Necessary to unpickle file!
from CutSharing.SDDP_utils import report_stats

if __name__ == '__main__':
    '''
    Implementation of Wasserstein uncertainty set based on taking
    the dual problem of the inner max problem that represents the
    worst-case expectation.
    '''
    load_algorithm_options()
    
    T, model_builder, random_builder, rnd_container_data, rnd_container_oos, r_dro, instance_name, instance_name_gen = load_hydro_data(
        'DUAL', 'DW')
    #options['cut_selector'] = SLACK_BASED_CUT_SELECTOR#SLACK_BASED_CUT_SELECTOR#LAST_CUTS_SELECTOR
    algo = SDDP(T,
                model_builder,
                random_builder,
                risk_measure=DistRobustWasserstein,
                norm=1,
                radius=0,
                data_random_container=rnd_container_data)
    lbs = algo.run(instance_name=instance_name)
    sim_result = algo.simulate_policy(rnd_container_oos)
    report_stats(sim_result.sims_ub)
    
    for r_DRO in dro_radii:
        instance_name = instance_name_gen(r_DRO)
        algo.change_dro_radius(r_DRO)
        lbs = algo.run(instance_name=instance_name)
        
        save_path = hydro_path + '/Output/DW_Dual/%s_LBS.pickle' % (instance_name)
        write_object_results(save_path, (algo.instance, lbs))
        
        sim_result = algo.simulate_policy(rnd_container_oos)
        save_path = hydro_path + '/Output/DW_Dual/%s_OOS.pickle' % (instance_name)
        write_object_results(save_path, sim_result)
        report_stats(sim_result.sims_ub)
    
    del (algo)
