from importFile import *
from pyomo.bilevel import *
from microgrid_Model import *
import os
import pandas as pd
import matplotlib.pyplot as plt
import optimizationModel,microgridStructure
from optimizationModel import *
'''Initialize a special case of microgrid'''
case = microgridStructure.case_IES
'''Load input data'''
microgrid_data = pd.read_excel('input_IES.xlsx')
'''Construct base model'''
optimalDispatch = DayAheadModel(microgrid_data,case,range(96))
for i in range(92):
    AddDayInSubModel(optimalDispatch, i, microgrid_data, case)
xfrm = TransformationFactory('bilevel.linear_mpec')
xfrm.apply_to(optimalDispatch)
xfrm = TransformationFactory('mpec.simple_disjunction')
xfrm.apply_to(optimalDispatch)
xfrm = TransformationFactory('gdp.bigm')
xfrm.apply_to(optimalDispatch, default_bigM=1000)
solver = SolverFactory('gurobi')
res = solver.solve(optimalDispatch)
case.update(optimalDispatch)
case.DCPowerFlow()

