'''
Created on Nov 17, 2017

@author: dduque
'''
import numpy as np
from gurobipy import *
from CutSharing import ZERO_TOL

class CutPool():
    
    def __init__(self, stage ):
        self.stage = stage
        self.pool = {}
        
        self._requires_ajustment = False
        
    
    def addCut(self, new_cut):
        if new_cut.recomputable_rhs:
            self._requires_ajustment = True
        self.pool[new_cut.name] = new_cut
    
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
        #Reference to the constraint
        self.ctrRef = m.addConstr( self.lhs >= self.rhs, self.name)
        
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
        
