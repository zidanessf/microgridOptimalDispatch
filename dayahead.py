from importFile import *
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
solver = SolverFactory('gurobi')
res = solver.solve(optimalDispatch,tee=True)
case.update(optimalDispatch)
case.DCPowerFlow()
df = pd.DataFrame()
T= optimalDispatch.T
for gt in case.getKey(microgridStructure.gasTurbine):
    df[gt] = [value(optimalDispatch.gt_power[gt,t]) for t in T]
for branch in case.graph.edges():
    nf = branch[0]
    nt = branch[1]
    df[str(nf) + ' to ' + str(nt) + ' power flow'] = case.graph.edge[nf][nt]['Power_Flow']
