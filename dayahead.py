from importFile import *
from pyomo.bilevel import *
from microgrid_Model import *
import os
import pandas as pd
import matplotlib.pyplot as plt
import optimizationModel,microgridStructure
from optimizationModel import *
def fix_master_var(mdl):
    mdl.sub.fix_master_var = ConstraintList()
    for (name, data) in mdl.component_map(active=True).items():
        if isinstance(data, Var):
            for v in data:
                mdl.sub.fix_master_var.add(expr=data[v]==0)
'''Initialize a special case of microgrid'''
case = microgridStructure.case_IES
'''Load input data'''
microgrid_data = pd.read_excel('input_IES.xlsx')
'''Construct base model'''
optimalDispatch = DayAheadModel(microgrid_data,case,range(96))
Tmpc = 3
subobj = 0
xfrm = TransformationFactory('bilevel.linear_mpec')
for i in range(Tmpc):
    AddDayInSubModel(optimalDispatch, i, microgrid_data, case)
    xfrm.apply_to(optimalDispatch.sub, options={'submodel': 'MPC_'+str(i)})
    temp = getattr(optimalDispatch.sub,'MPC_'+str(i))
    subobj += temp.obj_Cost(temp,0)
'''set subproblem objective'''
optimalDispatch.sub.o = Objective(expr=subobj)
'''fix master varibales'''
fix_master_var(optimalDispatch)
tempmodel = ConcreteModel()
xfrm = TransformationFactory('mpec.simple_disjunction')
xfrm.apply_to(optimalDispatch.sub)
xfrm = TransformationFactory('gdp.bigm')
xfrm.apply_to(optimalDispatch.sub, default_bigM=1000)
solver = SolverFactory('gurobi')
res = solver.solve(optimalDispatch.sub)

