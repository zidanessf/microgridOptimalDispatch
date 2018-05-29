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
optimalDispatch = optimizationModel.DayAheadModel(microgrid_data,case,range(32))
xfrm = TransformationFactory('bilevel.linear_mpec')
xfrm.apply_to(optimalDispatch)
xfrm = TransformationFactory('mpec.simple_disjunction')
xfrm.apply_to(optimalDispatch)
xfrm = TransformationFactory('gdp.bigm')
xfrm.apply_to(optimalDispatch, default_bigM=10000)
solver = SolverFactory('cplex')
res = solver.solve(optimalDispatch,tee=True)
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
