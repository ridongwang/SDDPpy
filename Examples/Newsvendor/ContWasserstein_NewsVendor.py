#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 16:13:29 2018

@author: dduque
"""

from gurobipy import *
import os
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
from CutSharing.RiskMeasures import DistRobustWassersteinCont
newsvendro_path = os.path.dirname(os.path.realpath(__file__))
from OutputAnalysis.SimulationAnalysis import SimResult, plot_sim_results
def solve_cont_wasserstain(N,K,L, xi_n, C, d, r, c, dro_radius):
    
    model = Model('ContWassersteinNewsVendor')
    model.params.OutputFlag = 0 
    model.params.Solver = 1
    x = model.addVar(lb=0,ub=GRB.INFINITY,obj=0,vtype=GRB.CONTINUOUS, name='x')
    #x = model.addVar(lb=3,ub=3,obj=0,vtype=GRB.CONTINUOUS, name='x')
    s = model.addVars(N,lb=-GRB.INFINITY,ub=GRB.INFINITY,obj=[1/len(N) for _ in range(n)], vtype=GRB.CONTINUOUS, name='s')
    lam = model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,obj=dro_radius,vtype=GRB.CONTINUOUS,name='lambda_var')
    gam = model.addVars(K,N,L,lb=0,ub=GRB.INFINITY,obj=0, vtype=GRB.CONTINUOUS, name='gamma_var')
    model.update()
    piece_ctrs = {}
    #Piece 0:
    for i in N:
        piece_ctrs[(0,i)] = model.addConstr(rhs=((c-r)*x+quicksum((d[l]-C[l]*xi_n[i])*gam[0,i,l] for l in L)), sense=GRB.GREATER_EQUAL, lhs=s[i], name='piece[0,%i]' %i)
    #Piece 1:
    for i in N:
        piece_ctrs[(1,i)] = model.addConstr(rhs=(c*x-r*xi_n[i]+quicksum((d[l]-C[l]*xi_n[i])*gam[1,i,l] for l in L)), sense=GRB.GREATER_EQUAL, lhs=s[i],name='piece[1,%i]' %i)
    
    norm_ctrs = {}
    for k in K:
        for i in N:
            linexp = LinExpr()
            if k == 0:
                linexp = quicksum(C[l]*gam[k,i,l] for l in L)
            else:
                linexp = quicksum(C[l]*gam[k,i,l] for l in L) + r
            norm_ctrs[(k,i,-1)] = model.addConstr(lhs=lam-linexp, sense=GRB.GREATER_EQUAL, rhs=0,name='InfNormNeg[%i,%i]' %(k,i)) #RESTA
            norm_ctrs[(k,i, 1)] =model.addConstr(lhs=lam+linexp, sense=GRB.GREATER_EQUAL, rhs=0,name='InfNormPos[%i,%i]' %(k,i)) #SUMA
    model.update()
    
    model.optimize()
    #===========================================================================
    # print(model.getObjective())
    # for c in model.getConstrs():
    #     print(c.ConstrName, model.getRow(c), c.Sense, c.RHS)
    #===========================================================================
    if model.status == GRB.OPTIMAL:
        print(model.ObjVal, x.X, [s[i].X for i in N])
        for c in model.getConstrs():
            if c.Pi  > 1E-8:
                print(c.ConstrName, c.Pi)
        print('orig supp' , xi_n)
        new_support = []
        pmf = []
        for (k,i) in piece_ctrs:
            if piece_ctrs[(k,i)].Pi > 1E-8:
                new_atom = xi_n[i]  + (norm_ctrs[(k,i,1)].Pi - norm_ctrs[(k,i,-1)].Pi)/piece_ctrs[(k,i)].Pi
                new_support.append(new_atom)
                pmf.append(piece_ctrs[(k,i)].Pi)
        new_support = np.array(new_support)
        pmf = np.array(pmf)
        supp_argsort = np.argsort(new_support) 
        pmf = pmf[supp_argsort]
        new_support.sort()

        for i in range(len(new_support)-1, 0, -1):
            if np.abs(new_support[i] -new_support[i-1])<1E-8:
                new_support = np.delete(new_support, obj=i)
                rep_prob = pmf[i]
                pmf =np.delete(pmf,obj =i)
                pmf[i-1] += rep_prob                
        print('new supp ', new_support)
        print('new pmf' , pmf)
        return x.X, new_support , pmf
    else:
        model.computeIIS()
        model.write("NewsVendorInf.ilp")
        os.system("open -a TextEdit filename Users/dduque/Desktop/NewsVendorInf.ilp")
        
def solve_cont_wasserstain2(N,K,L, xi_n, C, d, r, c, dro_radius):
    '''
    Different objective 
    '''
    model = Model('ContWassersteinNewsVendor')
    model.params.OutputFlag = 0 
    model.params.Solver = 1
    x = model.addVar(lb=0,ub=GRB.INFINITY,obj=0,vtype=GRB.CONTINUOUS, name='x')
    #x = model.addVar(lb=3,ub=3,obj=0,vtype=GRB.CONTINUOUS, name='x')
    s = model.addVars(N,lb=-GRB.INFINITY,ub=GRB.INFINITY,obj=[1/len(N) for _ in range(n)], vtype=GRB.CONTINUOUS, name='s')
    lam = model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,obj=dro_radius,vtype=GRB.CONTINUOUS,name='lambda_var')
    gam = model.addVars(K,N,L,lb=0,ub=GRB.INFINITY,obj=0, vtype=GRB.CONTINUOUS, name='gamma_var')
    model.update()
    piece_ctrs = {}
    #y=\max\left(\left(c-p_3\right)x,cx-\left(p_1D\ -40\right),cx-\left(p_2D\ -20\right),cx-p_3D\ \right)
    #Piece 0:
    for i in N:
        piece_ctrs[(0,i)] = model.addConstr(rhs=((c-r[0])*x+quicksum((d[l]-C[l]*xi_n[i])*gam[0,i,l] for l in L)), sense=GRB.GREATER_EQUAL, lhs=s[i], name='piece[0,%i]' %i)
    #Piece 1:
    for i in N:
        piece_ctrs[(1,i)] = model.addConstr(rhs=(c*x-r[0]*xi_n[i]+quicksum((d[l]-C[l]*xi_n[i])*gam[1,i,l] for l in L)), sense=GRB.GREATER_EQUAL, lhs=s[i],name='piece[1,%i]' %i)
    
    #Piece 2:
    for i in N:
        piece_ctrs[(2,i)] = model.addConstr(rhs=(c*x-r[1]*xi_n[i] +20 +quicksum((d[l]-C[l]*xi_n[i])*gam[2,i,l] for l in L)), sense=GRB.GREATER_EQUAL, lhs=s[i],name='piece[2,%i]' %i)
    
    #Piece 3:
    for i in N:
        piece_ctrs[(3,i)] = model.addConstr(rhs=(c*x-r[2]*xi_n[i] +40+quicksum((d[l]-C[l]*xi_n[i])*gam[3,i,l] for l in L)), sense=GRB.GREATER_EQUAL, lhs=s[i],name='piece[3,%i]' %i)
    
    
    norm_ctrs = {}
    for k in K:
        for i in N:
            linexp = LinExpr()
            if k == 0:
                linexp = quicksum(C[l]*gam[k,i,l] for l in L)
            else:
                linexp = quicksum(C[l]*gam[k,i,l] for l in L) + r[k-1]
            norm_ctrs[(k,i,-1)] = model.addConstr(lhs=lam-linexp, sense=GRB.GREATER_EQUAL, rhs=0,name='InfNormNeg[%i,%i]' %(k,i)) #RESTA
            norm_ctrs[(k,i, 1)] =model.addConstr(lhs=lam+linexp, sense=GRB.GREATER_EQUAL, rhs=0,name='InfNormPos[%i,%i]' %(k,i)) #SUMA
    model.update()
    
    model.optimize()
    #===========================================================================
    # print(model.getObjective())
    # for c in model.getConstrs():
    #     print(c.ConstrName, model.getRow(c), c.Sense, c.RHS)
    #===========================================================================
    if model.status == GRB.OPTIMAL:
        print(model.ObjVal, x.X, [s[i].X for i in N])
        for c in model.getConstrs():
            if c.Pi  > 1E-8:
                print(c.ConstrName, c.Pi)
        print('orig supp' , xi_n)
        new_support = []
        pmf = []
        for (k,i) in piece_ctrs:
            if piece_ctrs[(k,i)].Pi > 1E-8:
                new_atom = xi_n[i]  + (norm_ctrs[(k,i,1)].Pi - norm_ctrs[(k,i,-1)].Pi)/piece_ctrs[(k,i)].Pi
                new_support.append(new_atom)
                pmf.append(piece_ctrs[(k,i)].Pi)
        new_support = np.array(new_support)
        pmf = np.array(pmf)
        supp_argsort = np.argsort(new_support) 
        pmf = pmf[supp_argsort]
        new_support.sort()

        for i in range(len(new_support)-1, 0, -1):
            if np.abs(new_support[i] -new_support[i-1])<1E-8:
                new_support = np.delete(new_support, obj=i)
                rep_prob = pmf[i]
                pmf =np.delete(pmf,obj =i)
                pmf[i-1] += rep_prob                
        print('new supp ', new_support)
        print('new pmf' , pmf)
        return x.X, new_support , pmf
    else:
        model.computeIIS()
        model.write("NewsVendorInf.ilp")
        os.system("open -a TextEdit filename Users/dduque/Desktop/NewsVendorInf.ilp")
        
        
def create_obj(c,r):
    def news_vendor_obj(x, demand):
        return c*x -r*demand + r*np.maximum(0,demand-x)
    return news_vendor_obj

def create_obj2(c,r):
    #y=\max\left(\left(c-p_3\right)q,cq-\left(p_1x\ -40\right),cq-\left(p_2x\ -20\right),cq-p_3x\ \right)
    def news_vendor_obj(x, xi_n):
        return np.maximum((c-r[0])*x, np.maximum(c*x-r[0]*xi_n, np.maximum(c*x-r[1]*xi_n +20, c*x-r[2]*xi_n +40)))
    return news_vendor_obj


def test_out_of_sample(x_star, test_set, obj_fun, instance):
    objs = obj_fun(x_star,test_set)
    sim_res = SimResult(instance,objs)
    x_bar = np.mean(objs)
    std  = np.std(objs)
    if instance == None:
        print('%10s %10.4f %10.4f %10.4f' %('EV_Policy', x_star, x_bar, std))
    else:
        print('%10.4f %10.4f %10.4f %10.4f' %(instance['risk_measure_params']['radius'], x_star, x_bar, std))
    return sim_res

    
    
if __name__ == '__main__':
    np.random.seed(0)
    C = [-1,1]
    #r = 100
    r = [3,5,9]
    c = 2
    #quantile = (r-c)/r
    K = [0,1,2,3]
    L = [0,1]
    news_vendor_fun = create_obj(c,r)
    news_vendor_fun = create_obj2(c,r)
    #density = np.append(np.append(np.random.lognormal(0,1,size=1000) ,np.random.normal(10,1,size=2000)), np.random.normal(15,2,size=1000))
    density = np.abs(np.random.normal(50,25,size=2000))
    
    #===========================================================================
    # n, bins, patches = plt.hist(density, 100, normed=1, facecolor='green', alpha=0.75)
    # y = mlab.normpdf( bins, density.mean(), density.std())
    # l = plt.plot(bins, y, 'r--' , linewidth=1)
    # plt.show()
    #===========================================================================
    oos_sim = np.random.choice(density, size=1000, replace=False)
    #===========================================================================
    # n, bins, patches = plt.hist(oos_sim, 100, normed=1, facecolor='green', alpha=0.75)
    # y = mlab.normpdf( bins, oos_sim.mean(), oos_sim.std())
    # l = plt.plot(bins, y, 'r--' , linewidth=1)
    # plt.show()
    #===========================================================================
    xi_n = np.array([])
    for n in [10]:
        print('Solving n=',n)
        extra_poits = n- len(xi_n)

        #np.array([10,18,50,56])
        xi_n = np.append(xi_n , np.random.choice(density, size=extra_poits, replace=False))
        heights,bins = np.histogram(xi_n,bins=int(n/5))
        heights = heights/sum(heights)
        #plt.bar(xi_n,[1/n for _ in range(n)], color="blue", alpha=0.5)
        plt.bar(bins[:-1],heights,width=(max(bins) - min(bins))/len(bins), color="blue", alpha=0.8)
        xi_n.sort()
        #ev_policy = xi_n[int(quantile*n) - 1] 
        #oos_dem_bar = news_vendor_fun(ev_policy,oos_sim).mean()
        #test_out_of_sample(ev_policy, oos_sim, news_vendor_fun, None)
        N = [i for i in range(n)]
        d = [-xi_n.min() * 0.5, xi_n.max() * 1.5]
        
        sim_results = []
        #for dro_radius in [100]:
        for dro_radius in [10]:
            instance_name = 'NewVendor_N%i_CW' %(n)
            instance = {'risk_measure_params':{}}
            
            instance['risk_measure_params']['radius'] = dro_radius
            instance['risk_measure']= DistRobustWassersteinCont
            x_star, supp, pmf =   solve_cont_wasserstain2(N,K,L,xi_n, C, d, r,c, dro_radius)
            data_worst_case = np.random.choice(supp, size=1000000, p=pmf)
            heights,bins = np.histogram(data_worst_case,bins=int(n/5))
            heights = heights/sum(heights)
            plt.bar(bins[:-1],heights,width=(max(bins) - min(bins))/len(bins), color="red", alpha=0.5)
            #plt.bar(supp,pmf, color="red", alpha=0.5)
            #plt.show()
            sim_res  = test_out_of_sample(x_star,oos_sim,news_vendor_fun,instance)
            sim_results.append(sim_res)
        
        plot_sim_results(sim_results, newsvendro_path+'/Output/%s.pdf' %(instance_name), n, excel_file = False)
