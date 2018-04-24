from importFile import *
from pyomo.bilevel import *
from microgrid_Model import *
import os
import pandas as pd
import matplotlib.pyplot as plt
import optimizationModel,microgridStructure
'''Initialize a special case of microgrid'''
case = microgridStructure.case_PS
'''Load input data'''
microgrid_data = pd.read_excel('input_PS.xlsx')
'''Construct base model'''
optimalDispatch = optimizationModel.DayAheadModel(microgrid_data,case,range(96))
xfrm = TransformationFactory('bilevel.linear_mpec')
xfrm.apply_to(optimalDispatch)
xfrm = TransformationFactory('mpec.simple_disjunction')
xfrm.apply_to(optimalDispatch)
xfrm = TransformationFactory('gdp.bigm')
xfrm.apply_to(optimalDispatch, default_bigM=10000)
solvermngr = SolverManagerFactory('neos')
solver = SolverFactory('cplex')
res = solvermngr.solve(optimalDispatch,solver = solver)
# res = solver.solve(optimalDispatch,tee=True)
print(res)
case.update(optimalDispatch)
case.DCPowerFlow()
df = pd.DataFrame()
T= optimalDispatch.sub.T
for gt in case.getKey(microgridStructure.gasTurbine):
    df[gt] = [value(optimalDispatch.sub.gt_power[gt,t]) for t in T]
for wt in case.getKey(microgridStructure.PV):
    df[wt] = [value(optimalDispatch.wp[wt,t]) for t in T]
for branch in case.graph.edges():
    nf = branch[0]
    nt = branch[1]
    df[str(nf) + ' to ' + str(nt) + ' power flow'] = case.graph.edge[nf][nt]['Power_Flow']
df.to_excel('optimistic.xlsx')
# server = SolverManagerFactory('neos')
# opt = SolverFactory('cplex')
# server.solve(optimalDispatch,solver=opt)
# print('总运行成本：'+str(value(optimalDispatch.sub.obj_simple(optimalDispatch.sub))))
# '''Solve the base model'''
# solver = SolverFactory('bilevel_blp_global')
# options = {
#     'solver' : 'gurobi',
#     'bigM' : 1000
# }
# solver.solve(optimalDispatch)
# print('-------经济性最优----------')
# optimizationModel.retriveResult(microgrid_data,case,optimalDispatch)
# print('总运行成本：'+str(value(optimalDispatch.sub.obj_Economical(optimalDispatch.sub))))
# print('能效： '+str(1/value(optimalDispatch.sub.obj_Efficiency(optimalDispatch.sub))))