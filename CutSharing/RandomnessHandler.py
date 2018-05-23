'''
Created on Nov 19, 2017

@author: dduque
'''
import copy
import numpy as np
import scipy.sparse as sp

#NOT USED FOR THE MOMENT
class SamplePath():
    '''
    classdocs
    '''
    def __init__(self, sampler_func):
        '''
        Constructor
        '''
        self.sampledvalues = sampler_func
    
    def getStageRealizations(self, stage):
        assert stage in self.sampledvalues, 'Sample path has not being created'
        return self.sampledvalues[stage]

class RandomContainer:
    '''
    Class to store and handle all random elements
    Attributes:
        stage_vectors (list of StageRandomVector): A list 
        with a random vector per stage. Vector elements had
        specified dependencies with vectors of the previous
         stages if it applies.
    '''
    
    def __init__(self):
        self.stage_vectors = []
    
    def append(self, stageRndVector):
        '''
        Adds a random vector associated with a stage
        '''
        self.stage_vectors.append(stageRndVector)

    def getSamplePath(self,rnd_stream, ev = False):
        '''
        Generates a sample path. This is, a list of dictionaries
        where each dictionary corresponds to a different stage and 
        each dictionary has the realized value for each random element
        with the key being the name of the random element and the value
        if a float with the numeric value of the realization.
        
        Args:
            ev (bool): If the sample path is deterministic for expected value 
                analysis (Default value is false).
        
        Return:
            partial_sample_path (list of dic): a sample path including all stages.
        '''
        partial_sample_path = []
        if ev == True:
            partial_sample_path = [srv.get_ev_vector(partial_sample_path) for srv in self.stage_vectors]
        else:
            for srv in self.stage_vectors:
                stage_sample = srv.getSample(partial_sample_path,rnd_stream) 
                partial_sample_path.append(stage_sample)
        return partial_sample_path
    
    def getStageSample(self, t, partial_sample_path, rnd_stream):
        '''
        Generates a sample for a given stage and appends it to a partial path
        **** Modifies the partial_sample_path given as a parameter ****
        Args:
            t (int): stage number to be sample.
            partial_sample_path (list of dict): a sample path up to stage t-1.
        '''
        srv = self.stage_vectors[t]
        stage_sample = srv.getSample(partial_sample_path, rnd_stream) 
        partial_sample_path.append(stage_sample)
        
        
    def preprocess_randomness(self):
        for sv in self.stage_vectors:
            sv.preproces_randomness() 
    
    
    
    def enumerate_scenarios(self,):
        all_scenarios = []
        self.recursive_enum(all_scenarios, 0, [])
        return all_scenarios
        
        
    
    def recursive_enum(self, all_scenarios, stage, past):
        sv = self.stage_vectors[stage]
        for outcome in sv.outcomes:
            past.append(outcome)
            if stage < len(self.stage_vectors)-1:
                self.recursive_enum(all_scenarios, stage + 1, past)
            else:
                all_scenarios.append(copy.deepcopy(past))
            past.pop()
            
            
    def reset_to_nominal_dist(self):
        for srv in self.stage_vectors:
            srv.reset_to_nominal_dist()        
        
    
    #Aux methods to make the container iterable
    def __iter__(self):
        return (x for x in self.stage_vectors)
    def __getitem__(self, index):
        return self.stage_vectors[index]
    def __setitem__(self, key, value):
        self.stage_vectors[key] = value
    def __delitem__(self,key):
        self.stage_vectors.__delitem__(key)
    def __repr__(self):
        return [x for x in self.stage_vectors].__repr__()
    
class StageRandomVector:
    '''
    A class to represent the randomness of a stage
    Attributes:
        stage (int): Number of the corresponding stage.
        elements (dict (str,RandomElement)): Dictionary storing each random
            element of the vector by its name as a key.
        vector order (dict (str,int)): dictionary to map the names with a fixed
            order of the elements in the vector.
        outcomes (list of dict): a list of possible outcomes of the random vector.
            Each element of the list is a dictionary containing the numerical values
            for all the random elements of the vector. In the interstage dependent case,
            this values are the independent part only.  
        p (ndarray): vector with the probabilities of each outcome. This vector is modifiable
            depending on the risk measure.
        is_indipendent (bool): flag to distinguish between independent and dependent cases. 
    '''
    def __init__(self, stage):
        self.stage = stage
        self.elements = {}
        self.vector_order = {}
        self.outcomes = []
        self.ev_outcomes = []
        self.outcomes_dim = 0
        self.p = None
        self.p_copy = None
        self.is_independent = True 
        self._first_element_added = False
      
    def addRandomElememnt(self, ele_name, ele_outcomes , ele_prob=None):
        if self._first_element_added == False:
            self.outcomes_dim = len(ele_outcomes)
            for e in ele_outcomes:
                self.outcomes.append({})
                self.ev_outcomes.append({})   
            self._first_element_added = True
            if ele_prob==None:
                self.p = np.array([1/len(ele_outcomes) for x in ele_outcomes])
            else:
                self.p = np.array(ele_prob)
            self.p_copy = self.p.copy()
        else:
            assert len(self.outcomes)==len(ele_outcomes), "Random element with a different number of outcomes."
        
        self.vector_order[ele_name] = len(self.elements)#Defines order as it gets built. 
        self.elements[ele_name] = RandomElement(self.stage, ele_name , ele_outcomes)
        e_mean = np.mean(ele_outcomes)
        for (i,e) in enumerate(ele_outcomes):
            self.outcomes[i][ele_name] = e
            self.ev_outcomes[i][ele_name] = e_mean
        return self.elements[ele_name]
    
    def modifyOutcomesProbabilities(self, newp):
        assert len(newp)==len(self.outcomes)
        self.p = np.array(newp)
    
    def getSample(self, partial_sample_path, rnd_gen):
        '''
        Generates a sample for the associate stage random vector. 
        It receives the partial sample path in case the random vector
        is stagewise dependent.
        Attributes:
            partial_sample_path (list of dict): A list of the realizations
                of the previous stages. The dictionary has the same form as
                the method output.
        Return:
            satage_sample (dict of (str-float)): A dictionary containing the
                realized value of each random element.
        '''
        try:
            lucky_outcome = rnd_gen.choice([i for i in range(0,self.outcomes_dim)], p = self.p)
            stage_sample = {e.name:
                            e.getSample(lucky_outcome, e.get_depedencies_realization(partial_sample_path)) for e in self.elements.values()}
            return stage_sample
        except:
            print(self.p)
    
    def get_ev_vector(self,partial_sample_path):
        return {e.name: e.comput_element_ev(e.get_depedencies_realization(partial_sample_path)) for e in self.elements.values()}
    
    def getOutcomes(self, sample_path, ev):
        if self.is_independent:
            if ev ==True:
                return self.ev_outcomes
            return self.outcomes
        else:
            if ev ==True:
                raise 'Expected value mode not implemented for dependent case.'
            outcomes_copy = copy.deepcopy(self.outcomes)
            for e in self.elements.values():
                e_dependencies = e.get_depedencies_realization(sample_path)
                new_e_outcomes = e.compute_outcomes(e_dependencies)
                for i in range(0,len(new_e_outcomes)):
                    outcomes_copy[i][e.name] = new_e_outcomes[i]
                    
            return outcomes_copy
    
    def get_sorted_outcome(self, outcome_index):
        '''
        Retrieve a specific outcome and returned as a
        vector according the the vector_order mapping.
        Args:
            outcome_index (int): index of the outcome of interest
        Return:
            xi (ndarray): outcome in vector form
        '''
        xi = np.zeros(len(self.elements))
        for ele in self.vector_order:
            xi[self.vector_order[ele]] = self.outcomes[outcome_index][ele]
        return xi
        
    def preproces_randomness(self):
        vec_dim = len(self.elements)
        R_matrices = [] #Auto reg matrices
        created = False
        for e in self.elements.values():
            abs_order_e = self.vector_order[e.name]
            if e.has_dependencies():
                if created == False:
                    R_matrices = [sp.coo_matrix((vec_dim,vec_dim), dtype=np.float64) for de in e.dependencies]
                    created = True
                
                for (i, d_stage) in enumerate(e.dependencies):
                    #assert len(self.elements)==len(e.dependencies[d_stage])
                    for dep_e in e.dependencies[d_stage]:
                        abs_order_d = self.vector_order[dep_e]
                        #R_matrices[i][abs_order_e,abs_order_d] = e.dependencies[d_stage][dep_e]
                        R_matrices[i] =  R_matrices[i] + sp.coo_matrix(([e.dependencies[d_stage][dep_e]] , ([abs_order_e],[abs_order_d])), shape=(vec_dim,vec_dim ))
                 
                
                
                self.is_independent = False
                self.lag = e.lag
        self.autoreg_matrices = R_matrices       
                
    def reset_to_nominal_dist(self):
        self.p = self.p_copy.copy()
    
    def __repr__ (self):
        return 't=%i %s' %(self.stage, self.elements.keys().__repr__())
    
    
class RandomElement:
    '''
    Class to represent a random element
    Attributes:
        stage (int): Stage in which the random element realizes. 
        name (str): Name of the random element.
        outcomes (1D - ndarray): possible outcomes of the random element.
        current_sample_path(list of StageRandon)
        
        dependencies (dict of (int, dict of (str-float)): dependencies stored by 
            stage number (key of the outer dictionary) and the independent random
            element names (str key of the by random element). The value is the
            multiplicative coefficient.
    '''
    
    def __init__(self, stage ,name , rnd_outcomes):
        self.stage = stage
        self.name = name
        self.outcomes = np.array(rnd_outcomes)
        self.dependencyFunction = None
        self.dependencies = {}
        self.current_sample_path = None 
    def __repr__(self):
        return 't=%i>%s' %(self.stage, self.name)
    
    def has_dependencies(self):
        return True if len(self.dependencies)>0 else False
    
    def addDependecyFunction(self, dependency_coefs, depFunc):
        '''
        Specify a function for a dependent model.
        Args:
            dependency_coefs (dict of (int, dict of (str-float)): dependency coefficients
                stored by stage number (key of the outer dictionary) and the independent 
                random element names (str key of the by random element). The value 
                is the multiplicative coefficient.
            depFunc (::func::) a function that specifies how the dependency is computed. The
                    function receives to arguments:
                    rnd_element (RandomElement): the dependent random element.
                    realizations (dict): a dictionary of the realization with the same structure
                                         as the element dependencies dictionary.
                    The function returns a 1D-ndarry.
                    
        '''
        self.dependencyFunction = depFunc
        self.dependencies = dependency_coefs 
        min_stage_dep = min(x for x in self.dependencies)
        self.lag = self.stage - min_stage_dep
        
    def compute_outcomes(self, depedencies_realizations ):
        if self.dependencyFunction == None:
            return self.outcomes
        else:
            new_outcomes = self.dependencyFunction(self,depedencies_realizations)
            assert isinstance(new_outcomes, np.ndarray) and len(self.outcomes) == len(new_outcomes), \
                    "Dependency function returned outcomes that don't match the necessary dimensions (%i)" %(len(new_outcomes))
            return new_outcomes
    
    def comput_element_ev(self, depedencies_realizations):
        if self.dependencyFunction == None:
            return np.mean(self.outcomes)
        else:
            new_outcomes = self.dependencyFunction(self,depedencies_realizations)
            return np.mean(new_outcomes)
         
    def get_depedencies_realization(self, partial_sample_path):  
        '''
        Gets the necessary values for the element dependencies
        given a partial path up to stage = self.stage - 1.
        '''
        my_depedencies = {}
        for (stage, realization) in enumerate(partial_sample_path):
            if stage in self.dependencies:
                my_depedencies[stage] = {}
                for dname  in self.dependencies[stage]:
                    my_depedencies[stage][dname] = realization[dname]
        return my_depedencies           
                
                   
    def getSample(self, p, dependencies):
        '''
        Generates a sample of the random element give a probability vector
        and a list of dependencies.
        Attributes:
            p (1D ndarry): Array of probabilities for each outcome,
            
            dependencies (dict of dict): a dictionary of dependencies.
                The key of the outer dictionary is the stage corresponding
                to the dependency and the the value is another dictionary.
                The key of the inner dictionary is the name of the random
                element and the value is its numerical realization.
                e.g.: For a random element in stage 3
                dependencies={2:{'b_1':20, 'b_2':40},
                              1:{'b_1':5, 'b_2':3}}
                    So this random element depends on random elements
                    of stages 1 and 2.      
        '''
        generalized_outcomes = self.compute_outcomes(dependencies)
        return generalized_outcomes[p]
    
def AR1_depedency(rnd_ele , realizations):
    assert type(realizations) == dict, 'Realizations are not in the expected form.'
    assert len(realizations) == 1, "Dependency model has a lag > 1."
    AR1model = sum(realizations[k][j]*rnd_ele.dependencies[k][j] for k in realizations for j in realizations[k])
    return AR1model +  rnd_ele.outcomes

    
    
    
    
    
