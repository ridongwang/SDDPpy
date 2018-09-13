'''
Created on Sep 12, 2018

@author: dduque
'''
from CutSharing import options, LAST_CUTS_SELECTOR, load_algorithm_options
from CutSharing.SDDP_Alg import SDDP
from CutSharing.RiskMeasures import DistRobust, DiscreteWassersteinInnerSolver
from Utils.file_savers import write_object_results
from HydroModel import load_hydro_data, model_builder, random_builder, valley_chain_oos, random_builder_out_of_sample
from InstanceGen.ReservoirChainGen import read_instance, HydroRndInstance #Necessary to unpickle file!
    
if __name__ == '__main__':
    '''
    Implementation of Wasserstein uncertainty set based on Philpott et al.
    Inner worst-case expectation problem is solved in the backward pass to
    compute the cuts. This approach is regarded as Primal.
    '''
    load_algorithm_options()
    T, r_dro, instance_name = load_hydro_data()
    options['cut_selector'] = LAST_CUTS_SELECTOR
    algo = SDDP(T, model_builder, random_builder, risk_measure = DistRobust, dro_solver = DiscreteWassersteinInnerSolver,\
                dro_solver_params = {'norm': 1 , 'radius':r_dro})
    lbs = algo.run(instance_name=instance_name, dynamic_sampling=options['dynamic_sampling'])
    
    save_path = hydro_path+'/Output/DisceteWassersteinSingleCut/%s_LBS.pickle' %(instance_name)
    write_object_results(save_path, (algo.instance, lbs))      
    
    out_of_sample_rnd_cont = random_builder_out_of_sample(valley_chain_oos)
    sim_result = algo.simulate_policy(options['sim_iter'], out_of_sample_rnd_cont)
    save_path = hydro_path+'/Output/DisceteWassersteinSingleCut/%s_OOS.pickle' %(instance_name)
    write_object_results(save_path, sim_result)
          
    del(algo)
    
    
    
