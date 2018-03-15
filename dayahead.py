from importFile import *
from pyomo.bilevel import *
from microgrid_Model import *
import pandas as pd
import matplotlib.pyplot as plt
import optimizationModel,microgridStructure
'''Initialize a special case of microgrid'''
case = microgridStructure.MicrogridCase()
'''Load input data'''
microgrid_data = pd.read_excel('input.xlsx')
'''Construct base model'''
optimalDispatch = optimizationModel.DayAheadModel(microgrid_data,case,range(2))
'''Solve the base model'''
xfrm = TransformationFactory('bilevel.linear_mpec')
xfrm.apply_to(optimalDispatch)
optimalDispatch.pprint()
# solver = SolverFactory('bilevel_blp_global')
# solver.solve(optimalDispatch)
# print('-------经济性最优----------')
# print('总运行成本：'+str(value(optimalDispatch.sub.obj_Economical(optimalDispatch.sub))))
# print('能效： '+str(1/value(optimalDispatch.sub.obj_Efficiency(optimalDispatch.sub))))