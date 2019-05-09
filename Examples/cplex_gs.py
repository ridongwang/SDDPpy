'''
Created on Mar 6, 2019

@author: dduque
'''
import cplex as cpx


model  =cpx.Cplex()


x = list(model.variables.add([-1,-1], [0,0], [10,10], names=['x1','x2']))


model.linear_constraints.add([cpx.SparsePair([0,1],[1,1])], ['L'], [8], names=['c1'])
model.linear_constraints.add([cpx.SparsePair([0,1],[-2,0.5])], ['L'], [1], names=['c1'])


model.solve()

print(model.solution.get_float_qual