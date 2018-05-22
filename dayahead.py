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
microgrid_data = pd.read_excel('input_IES.xlsx')

'''Construct base model'''
T = 24
Lmpc = 2
Tmpc = T - Lmpc - 1
OD = DayAheadModel(microgrid_data,case,range(T))
# get the KKT conditions
subobj = 0
fix_all_var(OD)
for i in range(Tmpc):
    AddDayInSubModel(OD, i, microgrid_data, case,range(Lmpc))
    temp = getattr(OD.sub, 'MPC_' + str(i))
    xfrm = TransformationFactory('bilevel.linear_mpec')
    if i >= 1:
        prev = getattr(OD.sub, 'MPC_' + str(i-1))
        fix_all_var(prev)
        xfrm.apply_to(OD.sub,options={'submodel':'MPC_' + str(i)})
        unfix_all_vars(prev)
    else:
        xfrm.apply_to(OD.sub,options={'submodel':'MPC_' + str(i)})
    if i == Tmpc - 1:
        subobj += temp.obj_Cost(temp, range(Lmpc))
    else:
        subobj += temp.obj_Cost(temp, 0)
    print('MPC_{0} created'.format(i))
# set sub-problem objective///fix master vars///transform the sub-problem
OD.sub.o = Objective(expr=subobj,sense=maximize)
OD.sub.obj_copy = subobj
transform_sub(OD.sub)
unfix_all_vars(OD)
transform_master(OD)
#'''THIS IS A TEST'''
# #Transformation
# transform_master(OD)
# '''the base model is constructed'''
# OD.sub.deactivate()
# solver = SolverFactory('gurobi')
# res = solver.solve(OD)
# '''The KKT&G Algorithm begins'''
# OD.sub.activate()
# fix_all_var(OD)
# solver = SolverFactory('gurobi')
# res = solver.solve(OD.sub, tee=True,  # stream the solver output
#                    keepfiles=True,  # print the MILP file for examination
#                    symbolic_solver_labels=True)  # fix变量之后，submodel可以直接求解
lb = - np.inf
ub = np.inf
NumIter = 1
while 1:
    # solve the master problem
    print('Iteration num {0}'.format(NumIter))
    solver = SolverFactory('cplex')
    OD.sub.deactivate()
    if NumIter >= 2:
        del OD.objective
        OD.objective = Objective(rule=lambda mdl: mdl.obj_Economical(mdl) + mdl.eta)
    res = solver.solve(OD,tee=True)
    # print([value(OD.utility_power[t]) for t in OD.T])
    if NumIter >= 2:
        lb = value(OD.objective)
        print('各子问题eta')
        for i in range(1,NumIter):
            try:
                sub_p = getattr(OD,'KKTG'+str(i))
            except Exception:
                continue
            print('子问题{0}'.format(i))
            temp_obj = value(sub_p.obj_copy)
            print(temp_obj)
            print([value(sub_p.pv[t]) for t in OD.T])
            if i == 1:
                pv_ref = [value(sub_p.pv[t]) for t in OD.T]
            else:
                if sum(abs((value(sub_p.pv[t])-pv_ref[t])) for t in OD.T) <= 0.1:
                    OD.del_component('KKTG'+str(i))
                    OD.del_component('opt_cut_' + str(i))
        del pv_ref
    # if res.solver.termination_condition == TerminationCondition.optimal:
    #     lb = value(OD.objective)
    # elif res.solver.termination_condition == TerminationCondition.unbounded:
    #     lb = - np.inf
    if NumIter == 1:
        print('master problem optimal value is {0}'.format(lb))
        print('master problem optimal eta is None')
        print('the lower bound is updated to {0}'.format(-np.inf))
    else:
        print('master problem optimal value is {0}'.format(lb - value(OD.eta)))
        print('master problem optimal eta is {0}'.format(value(OD.eta)))
        print('the lower bound is updated to {0}'.format(lb))
    # solve the sub problem
    OD.sub.activate()
    fix_all_var(OD) #注意查看文档，fix变量的方法,fix master variables
    solver = SolverFactory('cplex')
    res = solver.solve(OD.sub,tee=True#stream the solver output
                         #print the MILP file for examination
                        ) #fix变量之后，submodel可以直接求解
    # res = solver.solve(OD.sub)
    if res.solver.termination_condition == TerminationCondition.optimal:
        ub = min(ub,value(OD.obj_Economical(OD)) + value(OD.sub.o))
    elif res.solver.termination_condition == TerminationCondition.unbounded:
        ub = ub
    # add kkt cuts to master problem
    print('the sub-obj is {0}'.format(value(OD.sub.o)))
    print('the upper bound is updated to {0}'.format(ub))
    # print([value(OD.sub.o)])
    # print([value(OD.sub.MPC_0.det_cs_cold['CS_1',t]) for t in range(4)])
    print('-----------------------')
    print('THE GAP IS {0}%'.format(100*(ub-lb)/(ub+lb)))
    print('-----------------------')
    if (ub-lb)/(ub+lb) <= 0.005:
        print('CONVERGED!')
        break
    add_kkt_cuts(OD, NumIter)
    unfix_all_vars(OD)
    NumIter += 1
