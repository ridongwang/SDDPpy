'''
Created on Jun 13, 2018

@author: dduque
'''



def print_model(model):
    '''
    Print out a gurobi model
    '''
    for c in model.getConstrs():
        print(c.ConstrName, model.getRow(c), c.Sense, c.RHS)