from importFile import *
from pyomo.bilevel import *
from microgrid_Model import *
import os
import pandas as pd
import matplotlib.pyplot as plt
from copy import deepcopy
import optimizationModel,microgridStructure
from optimizationModel import *
def fix_master_var(mdl,refmdl):
    mdl.fix_master_var = ConstraintList()
    for (name, data) in mdl.component_map(active=True).items():
        if isinstance(data, Var):
            for v in data:
                ref_v = getattr(refmdl,name)
                try:
                    mdl.fix_master_var.add(expr=data[v]==value(ref_v[v])) #固定主问题变量
                except Exception:
                    pass
'''Initialize a special case of microgrid'''
case = microgridStructure.case_IES
'''Load input data'''
microgrid_data = pd.read_excel('input_IES.xlsx')
'''Construct base model'''
optimalDispatch = DayAheadModel(microgrid_data,case,range(96))
subproblem = ConsFreeModel(microgrid_data,case,range(96))
Tmpc = 3
subobj = 0
xfrm = TransformationFactory('bilevel.linear_mpec')
for i in range(Tmpc):
    AddDayInSubModel(optimalDispatch, i, microgrid_data, case)
    xfrm.apply_to(optimalDispatch.sub, options={'submodel': 'MPC_'+str(i)})
    temp = getattr(optimalDispatch.sub, 'MPC_' + str(i))
    temp.deactivate()
for i in range(Tmpc):
    AddDayInSubModel(subproblem, i, microgrid_data, case)
    xfrm.apply_to(subproblem.sub, options={'submodel': 'MPC_'+str(i)})
    temp = getattr(subproblem.sub,'MPC_'+str(i))
    subobj += temp.obj_Cost(temp,0)
    temp.deactivate()
'''set subproblem objective'''
subproblem.o = Objective(expr=subobj)
'''solve the master problem'''
optimalDispatch.sub.deactivate()
xfrm = TransformationFactory('mpec.simple_disjunction')
xfrm.apply_to(optimalDispatch)
xfrm = TransformationFactory('gdp.bigm')
xfrm.apply_to(optimalDispatch, default_bigM=1000)
solver = SolverFactory('gurobi')
res = solver.solve(optimalDispatch)
'''fix master varibales'''
# fix_master_var(subproblem,optimalDispatch)
xfrm = TransformationFactory('mpec.simple_disjunction')
xfrm.apply_to(subproblem.sub)
xfrm = TransformationFactory('gdp.bigm')
xfrm.apply_to(subproblem.sub, default_bigM=1000)
solver = SolverFactory('gurobi')
res = solver.solve(subproblem)
