'''
Created on Nov 19, 2017

@author: dduque
'''
import copy
import numpy as np
np.random.seed(0)

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

    def getSamplePath(self):
        '''
        Generates a sample path. This is, a list of dictionaries
        where each dictionary corresponds to a different stage and 
        each dictionary has the realized value for each random element
        with the key being the name of the random element and the value
        if a float with the numeric value of the realization.
        '''
        partial_sample_path = []
        for srv in self.stage_vectors:
            stage_sample = srv.getSample(partial_sample_path) 
            partial_sample_path.append(stage_sample)
        return partial_sample_path
    
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
    A class to represent the randomnes of a stage
    Attributes:
        stage (int): Number of the corresponding stage.
    '''
    def __init__(self, stage):
        self.stage = stage
        self.elements = {}
        self.vector_order = {}
        self.outcomes = []
        self.outcomes_dim = 0
        self._first_element_added = False
        self.p = None
        self.is_independent = True 
      
    def addRandomElememnt(self, ele_name, ele_outcomes , ele_prob=None):
        if self._first_element_added == False:
            self.outcomes_dim = len(ele_outcomes)
            for e in ele_outcomes:
                self.outcomes.append({})   
            self._first_element_added = True
            if ele_prob==None:
                self.p = np.array([1/len(ele_outcomes) for x in ele_outcomes])
            else:
                self.p = np.array(ele_prob)
        else:
            assert len(self.outcomes)==len(ele_outcomes), "Random element with a different number of outcomes."
        
        self.vector_order[ele_name] = len(self.elements)#Defines order as it gets built. 
        self.elements[ele_name] = RandomElement(self.stage, ele_name , ele_outcomes)
        for (i,e) in enumerate(ele_outcomes):
            self.outcomes[i][ele_name] = e
        return self.elements[ele_name]
    
    def modifyOutcomesProbabilities(self, newp):
        assert len(newp)==len(self.outcomes)
        self.p = np.array(newp)
    
    def getSample(self, partial_sample_path):
        '''
        Generates a sample for the associate stage random vector. 
        It receives the partial sample path in case the random vector
        is stagewise dependent.
        Attributes:
            partial_sample_path (list of dict): A list of the realizations
                of the previous stages. The dictonary has the same form as
                the method output.
        Return:
            satage_sample (dict of (str-float)): A dictionary containing the
                realized value of each random element.
        '''
        lucky_outcome = np.random.choice([i for i in range(0,self.outcomes_dim)], p = self.p)
        stage_sample = {e.name:
                        e.getSample(lucky_outcome, e.get_depedencies_realization(partial_sample_path)) for e in self.elements.values()}
        return stage_sample
    
    def getOutcomes(self, sample_path):
        if self.is_independent:
            return self.outcomes
        else:
            outcomes_copy = copy.deepcopy(self.outcomes)
            for e in self.elements.values():
                e_dependencies = e.get_depedencies_realization(sample_path)
                new_e_outcomes = e.compute_outcomes(e_dependencies)
                for i in range(0,len(new_e_outcomes)):
                    outcomes_copy[i][e.name] = new_e_outcomes[i]
                    
            return outcomes_copy
    
    
    def preproces_randomness(self):
        vec_dim = len(self.elements)
        R_matrices = [] #Auto reg matrices
        created = False
        for e in self.elements.values():
            abs_order_e = self.vector_order[e.name]
            if e.has_dependencies():
                if created == False:
                    R_matrices = [np.zeros((vec_dim,vec_dim)) for de in e.dependencies]
                    created = True
                
                for (i, d_stage) in enumerate(e.dependencies):
                    assert len(self.elements)==len(e.dependencies[d_stage])
                    for dep_e in e.dependencies[d_stage]:
                        abs_order_d = self.vector_order[dep_e]
                        R_matrices[i][abs_order_e,abs_order_d] = e.dependencies[d_stage][dep_e]
                 
                
                
                self.is_independent = False
                self.lag = e.lag
        self.autoreg_matrices = R_matrices       
                
        
    
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
        #return np.random.choice(generalized_outcomes, p=p)
    
def AR1_depedency(rnd_ele , realizations):
    assert type(realizations) == dict, 'Realizations are not in the expected form.'
    assert len(realizations) == 1, "Dependency model has a lag > 1."
    AR1model = sum(realizations[k][j]*rnd_ele.dependencies[k][j] for k in realizations for j in realizations[k])
    return AR1model +  rnd_ele.outcomes

    
    
    
    
    