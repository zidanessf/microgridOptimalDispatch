from importFile import *
from pyomo.bilevel import *
from microgrid_Model import *
import os
import pandas as pd
import matplotlib.pyplot as plt
import optimizationModel,microgridStructure
'''Initialize a special case of microgrid'''
case = microgridStructure.MicrogridCase()
'''Load input data'''
microgrid_data = pd.read_excel('input.xlsx')
'''Construct base model'''
optimalDispatch = optimizationModel.DayAheadModel(microgrid_data,case,range(2))
xfrm = TransformationFactory('bilevel.linear_mpec')
xfrm.apply_to(optimalDispatch)
xfrm = TransformationFactory('mpec.simple_disjunction')
xfrm.apply_to(optimalDispatch)
xfrm = TransformationFactory('gdp.bigm')
xfrm.apply_to(optimalDispatch, default_bigM=1000)
solver = SolverFactory('gurobi')
solver.solve(optimalDispatch)
print('总运行成本：'+str(value(optimalDispatch.sub.obj_Economical(optimalDispatch.sub))))
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