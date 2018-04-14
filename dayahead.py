from importFile import *
from pyomo.bilevel import *
from microgrid_Model import *
import pandas as pd
import matplotlib.pyplot as plt
import microgridStructure
import numpy as np
from optimizationModel import *

'''Initialize a special case of microgrid'''
case = microgridStructure.case_IES
'''Load input data'''
microgrid_data = pd.read_excel('IES_SIMPLE.xlsx')

'''Construct base model'''
OD = DayAheadModel(microgrid_data,case,range(5))
Tmpc = 2
Lmpc = 4
xfrm = TransformationFactory('bilevel.linear_mpec')
subobj = 0
# get the KKT conditions
for i in range(Tmpc):
    AddDayInSubModel(OD, i, microgrid_data, case)
    xfrm.apply_to(OD.sub, options={'submodel': 'MPC_'+str(i)})
    temp = getattr(OD.sub, 'MPC_' + str(i))
    subobj += temp.obj_Cost(temp,0)
# set sub-problem objective
OD.sub.o = Objective(expr=subobj,sense=maximize)
#Transformation
xfrm = TransformationFactory('mpec.simple_disjunction')
xfrm.apply_to(OD.sub)
xfrm = TransformationFactory('gdp.bigm')
xfrm.apply_to(OD.sub, default_bigM=100000) #Big-M过小会导致不可行
xfrm = TransformationFactory('mpec.simple_disjunction')
xfrm.apply_to(OD)
xfrm = TransformationFactory('gdp.bigm')
xfrm.apply_to(OD, default_bigM=1000)
'''the base model is constructed'''

'''The KKT&G Algorithm begins'''
lb = - np.inf
ub = np.inf
NumIter = 1
while 1:
    # solve the master problem
    print('Iteration num {0}'.format(NumIter))
    solver = SolverFactory('gurobi')
    OD.sub.deactivate()
    res = solver.solve(OD)
    if res.solver.termination_condition == TerminationCondition.optimal:
        lb = value(OD.objective)
    elif res.solver.termination_condition == TerminationCondition.unbounded:
        lb = - np.inf
    print('master problem optimal value is {0}'.format(lb - value(OD.eta)))
    print('master problem optimal eta is {0}'.format(value(OD.eta)))
    print('the lower bound is updated to {0}'.format(lb))
    print('the fuck heat_ex is ')
    print([value(OD.heat_ex[t]) for t in range(5)])
    # solve the sub problem
    OD.sub.activate()
    fix_all_var(OD) #注意查看文档，fix变量的方法,fix master variables
    solver = SolverFactory('gurobi')
    res = solver.solve(OD.sub,tee=True, #stream the solver output
                        keepfiles=True, #print the MILP file for examination
                        symbolic_solver_labels=True) #fix变量之后，submodel可以直接求解
    print(res)
    if res.solver.termination_condition == TerminationCondition.optimal:
        ub = min(ub,value(OD.obj_Economical(OD)) + value(OD.sub.o))
    elif res.solver.termination_condition == TerminationCondition.unbounded:
        ub = ub
    # add kkt cuts to master problem
    print('the upper bound is updated to {0}'.format(ub))
    add_kkt_cuts(OD, NumIter)
    unfix_all_vars(OD)
    NumIter += 1
