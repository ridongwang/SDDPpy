'''
Created on Jun 13, 2018

@author: dduque
'''



def print_model(m):
    '''
    Print out a gurobi model
    '''
    print(m.getConstrs())
    for c in m.getConstrs(): print(c.ConstrName, m.getRow(c) , '  ', c.Sense, '  ', c.RHS)
    for v in m.getVars(): print(v.varname, ' '  , v.lb , '  ---  ', v.ub)