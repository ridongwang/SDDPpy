'''
Created on Sep 13, 2018

@author: dduque
'''
from CutSharing import options, LAST_CUTS_SELECTOR, load_algorithm_options,\
    SLACK_BASED_CUT_SELECTOR
from CutSharing.SDDP_Alg import SDDP
from CutSharing.RiskMeasures import DistRobustWasserstein
from Utils.file_savers import write_object_results
from HydroModel import load_hydro_data, model_builder, random_builder, valley_chain_oos, random_builder_out_of_sample
from InstanceGen.ReservoirChainGen import read_instance, HydroRndInstance #Necessary to unpickle file!
from __init__ import hydro_path

if __name__ == '__main__':
    '''
    Implementation of Wasserstein uncertainty set based on taking
    the dual problem of the inner max problem that represents the
    worst-case expectation.
    '''
    load_algorithm_options()
    T, r_dro, instance_name, out_of_sample_rnd_cont = load_hydro_data('DUAL', 'DW')
    options['cut_selector'] = LAST_CUTS_SELECTOR#SLACK_BASED_CUT_SELECTOR#LAST_CUTS_SELECTOR
    algo = SDDP(T, model_builder, random_builder, risk_measure = DistRobustWasserstein , norm = 1 , radius = r_dro)
    lbs = algo.run(instance_name=instance_name, dynamic_sampling=options['dynamic_sampling'])                                                              
    
    save_path = hydro_path+'/Output/DW_Dual/%s_LBS.pickle' %(instance_name)
    write_object_results(save_path, (algo.instance, lbs))
    
    sim_result = algo.simulate_policy(options['sim_iter'], out_of_sample_rnd_cont)
    save_path = hydro_path+'/Output/DW_Dual/%s_OOS.pickle' %(instance_name)
    write_object_results(save_path, sim_result)
    
    del(algo)
    
    
