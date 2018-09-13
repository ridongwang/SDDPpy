'''
Created on Nov 17, 2017

@author: dduque
'''
import numpy as np
from gurobipy import *
from CutSharing import ZERO_TOL, options, LAST_CUTS_SELECTOR, SLACK_BASED_CUT_SELECTOR
from abc import ABC, abstractmethod

class CutPool():
    
    def __init__(self, stage ):
        self.stage = stage
        self.pool = {}
        self.pool_order = []
        self.cut_selector = None
        self._requires_ajustment = False
        
        if options['cut_selector'] == LAST_CUTS_SELECTOR:
            self.cut_selector = LastCutsSelector()
        elif options['cut_selector'] == SLACK_BASED_CUT_SELECTOR:
            self.cut_selector = SlackBasedCutSelector()
        else:
            self.cut_selector = None
            
        
        
    
    def addCut1(self, model, new_cut):
        if new_cut.recomputable_rhs:
            self._requires_ajustment = True
        assert new_cut.name not in self.pool
        self.pool[new_cut.name] = new_cut
        self.pool_order.append(new_cut.name)
        if self.cut_selector != None:
            self.cut_selector.select_cuts(model, self.pool, self.pool_order)
    
    def addCuts(self, model, new_cuts):
        '''
        Add cut to the pool manager
        
        Params:
            model(GRBModel): Model that contain the cuts
            new_cuts (list of cut): list with the cuts to be added.
        '''
        for new_cut in new_cuts:
            if new_cut.recomputable_rhs:
                self._requires_ajustment = True
            assert new_cut.name not in self.pool
            self.pool[new_cut.name] = new_cut
            self.pool_order.append(new_cut.name)

        if self.cut_selector != None:
            self.cut_selector.add_recent_cuts([(c.name, c.is_active) for c in new_cuts])
            self.cut_selector.select_cuts(model, self.pool, self.pool_order)
    
    
    def needs_update(self):
        return self._requires_ajustment
    
    
    def get_non_zero_duals(self):
        tup_ind = []
        duals = []
        for ctr_name in self.pool:
            opt_cut = self.pool[ctr_name]
            if np.abs(opt_cut.ctrRef.Pi) > ZERO_TOL:
                tup_ind.append((opt_cut.cut_id, opt_cut.outcome))
                duals.append(opt_cut.ctrRef.Pi)
        return tup_ind, duals
                
    
    def __len__(self):
        return len(self.pool)
    def __iter__(self):
        return (x for x in self.pool.values())
        
class Cut():
    '''
    Structure for a cut
    '''
    def __init__(self,sp, var_coeffs, intercept, cut_id , stagewise_ind = True, ind_rhs = None, dep_rhs_vector = None, outcome = 0):
        '''
        Args:
            var_coeffs (dict of reals): coefficient of each variable involved in the cut where variable name is the key.
            vars (dict of GRBVar): dictionary of decision variables where the variable name is the key.
            
        '''
        m = sp.model
        stage = sp.stage
        self.name = 'cut[%i,%i,%i]' %(stage,cut_id,outcome)
        self.lhs = quicksum(-var_coeffs[vn]*m.getVarByName(vn) for vn in var_coeffs) 
        self.lhs.add(sp.oracle[outcome])
        self.recomputable_rhs = (stagewise_ind==False)
        self.rhs = intercept
        self.ind_rhs = ind_rhs
        self.dep_rhs_vector = dep_rhs_vector
        #=======================================================================
        # self.ctrRef = m.addConstr( self.lhs >= self.rhs, self.name)
        # self.is_active  = True
        #=======================================================================
        #Reference to the constraint
        if self.lhs.getValue() > self.rhs - ZERO_TOL:
            self.is_active  = False
        else:
            #print('New cuts: ' , self.name)
            self.is_active  = True
            #Reference to the constraint
            self.ctrRef = m.addConstr( self.lhs >= self.rhs, self.name)
            #print(self.lhs.getValue()  , self.rhs)
        
        #==============================#
        # Extra information for dual retrieval
        self.cut_id = cut_id
        self.outcome = outcome

        

    def adjust_intercept(self, omega_last):
        '''
        omega_last (1D ndarray): current ancestor scenario.
        '''
        dep_rhs = self.dep_rhs_vector.dot(omega_last).squeeze()
        new_rhs = dep_rhs + self.ind_rhs 
        self.ctrRef.RHS = new_rhs
        

class CutSelector(ABC):
    @abstractmethod
    def __init__(self):
        self.active = []
        self.unactive = []
    
    @abstractmethod
    def select_cuts(self, model, pool, pool_order):
        '''
        Makes a selection of the cuts to consider in the model.
        Updates local active and unactive list as well as the 
        status of the cut (cut.is_active and cut.ctrRef)
        '''
        raise 'Unimplmented method in a cut selection class.'
    
    def enforce_all(self,model,pool):
        for u in self.unactive:
            pool[u].ctrRef = model.addConstr( pool[u].lhs  >= pool[u].rhs, pool[u].name)
            self.active.append(u)
        self.unactive.clear()
        
    def add_recent_cuts(self,c_info):
        for (c_name,c_is_active) in c_info:
            if c_is_active:
                self.active.append(c_name)
            else:
                self.unactive.append(c_name)

class LastCutsSelector(CutSelector):
    '''
    Cut selector based on the last added cuts. The number of cuts
    is set in the algorithm options.
    '''
    def __init__(self):
        super().__init__()
        
    def select_cuts(self, model, pool, pool_order):
        N_max = options['max_cuts_last_cuts_selector']
        if len(self.active)<= N_max or len(self.active)==0:
            pass
        else:
            iters = 0
            while len(self.active)> N_max:
                a = self.active[0]
                if pool[a].lhs.getValue() > pool[a].rhs + 1E-6: #todo: use other parameter
                    pool[a].is_active = False
                    model.remove(pool[a].ctrRef)
                    pool[a].ctrRef = None
                    self.unactive.append(a)
                    self.active.pop(0)
                else:
                    iters +=1
                    self.active.append(self.active.pop(0))
                if iters >= len(self.active):
                    '''
                    Cardinality of selected cuts can't be satisfied.
                    Augmenting the maximum number of cuts
                    '''
                    options['max_cuts_last_cuts_selector'] = int(1.5*options['max_cuts_last_cuts_selector'])
                    print('max_cuts_last_cuts_selector -- > ' , options['max_cuts_last_cuts_selector'])
                    break
                
                
        

class  SlackBasedCutSelector(CutSelector):
    '''
    Cut selector based on binding cuts. At each iteration
    statistics on the cuts are updated to keep track of the
    number of times a cut is non-binding and is removed from
    the problem after a threshold is exceeded. 
    
    Attributes:
        active_stats (dict of (str,int)): count the number of times a cut 
            has been non-binding.
    '''
    def __init__(self):
        super().__init__()
        self.active_stats = {}
        self.track_ini=0
        self.track_end=0

    def select_cuts(self, model, pool, pool_order):
        
        #Update stats 
        tracked = int(len(self.active))
        ini_changed = False
        for i in range(self.track_ini, tracked):
            a = self.active[i]
            if (a in self.active_stats) == False:
                self.active_stats[a]=0
            else:
                #Count unactive times
                if pool[a].lhs.getValue() >= pool[a].rhs + options['slack_cut_selector']:
                    self.active_stats[a] +=1
                    if self.active_stats[a] >= options['slack_num_iters_cut_selector'] and ini_changed==False:
                        self.track_ini= i
                    else:
                        ini_changed=True
                else:    
                    self.active_stats[a] = 0 
                
            
        N_max = options['max_cuts_slack_based']
        if len(self.active)<= N_max or len(self.active)==0:
            pass #Not cut management requiered
        else:
            #Check active cuts (removing cuts from subproblems)
            hay_cut = len(self.active)
            new_active = []
            for i in range(tracked):
                a = self.active[i]
                if self.active_stats[a] >= options['slack_num_iters_cut_selector']:
                    pool[a].is_active = False
                    model.remove(pool[a].ctrRef)
                    pool[a].ctrRef = None
                    self.unactive.append(a)
                else:
                    new_active.append(a)
            new_active.extend(self.active[tracked:])
            self.active = new_active
            
            #Check unactive cuts (adding cuts)
            new_unactive = []
            for u in self.unactive:
                if pool[u].lhs.getValue()  < pool[u].rhs:
                    self.active.append(u)
                    pool[u].ctrRef = model.addConstr( pool[u].lhs  >= pool[u].rhs, pool[u].name)
                    pool[u].is_active = True
                    self.active_stats[u] = 0
                else: 
                    new_unactive.append(u)
            self.unactive = new_unactive
            self.track_ini=0
            print(self.active[0], hay_cut, len(self.active))
            if len(self.active)>N_max:
                options['max_cuts_slack_based'] = int(1.5*options['max_cuts_slack_based'])
                print('max_cuts_slack_based -- > ' , options['max_cuts_slack_based'])
            