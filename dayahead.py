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
optimalDispatch = optimizationModel.DayAheadModel(microgrid_data,case,range(96))
'''Solve the base model'''
xfrm = TransformationFactory('bilevel.linear_mpec')
xfrm.apply_to(optimalDispatch)
xfrm = TransformationFactory('bilevel.linear_mpec')
xfrm.apply_to(optimalDispatch)
solver = SolverFactory('gurobi')
solver.solve(optimalDispatch)
print('-------经济性最优----------')
print('总运行成本：'+str(value(optimalDispatch.sub.obj_Economical(optimalDispatch))))
print('能效： '+str(1/value(optimalDispatch.sub.obj_Efficiency(optimalDispatch))))
raise Exception
c11 = value(optimalDispatch.obj_Economical(optimalDispatch))
c12 = value(optimalDispatch.obj_Efficiency(optimalDispatch))
del optimalDispatch.objective
optimalDispatch.objective = Objective(rule=optimalDispatch.obj_Efficiency)
xfrm = TransformationFactory('gdp.chull')
xfrm.apply_to(optimalDispatch)
solver = SolverFactory('glpk')
solver.solve(optimalDispatch)
print('-------能效最优----------')
print('总运行成本：'+str(value(optimalDispatch.obj_Economical(optimalDispatch))))
print('能效： '+str(1/value(optimalDispatch.obj_Efficiency(optimalDispatch))))
c21 = value(optimalDispatch.obj_Economical(optimalDispatch))
c22 = value(optimalDispatch.obj_Efficiency(optimalDispatch))
del optimalDispatch.objective
k1 = c21 - c11
k2 = - c22 + c12
optimalDispatch.objective = Objective(rule= lambda mdl:k2*optimalDispatch.obj_Economical(mdl)+k1*optimalDispatch.obj_Efficiency(mdl))
xfrm = TransformationFactory('gdp.chull')
xfrm.apply_to(optimalDispatch)
solver = SolverFactory('glpk')
solver.solve(optimalDispatch)
print('-------' + str(k2)+'*经济性+' + str(k1) +'*能效最优----------')
print('总运行成本：'+str(value(optimalDispatch.obj_Economical(optimalDispatch))))
print('能效： '+str(1/value(optimalDispatch.obj_Efficiency(optimalDispatch))))
raise Exception('优化结束')
'''Retrieve the result'''
result = optimizationModel.retriveResult(microgrid_data,case,optimalDispatch)
result.to_excel('output.xlsx')
optimizationModel.extendedResult(result)
'''Interaction process'''
'''Electric DR'''
writer = pd.ExcelWriter('DayAhead Electric DR.xlsx')
peak = range(72, 76)
(max_model,max_amount) = optimizationModel.getMaxAmount(optimalDispatch,case,peak=peak,amount = [3000]*len(peak),mode='E')
for e in [1,0.8,0.6,0.4,0.2]:
    amount = [e * max_amount[t] for t in range(max_amount.__len__())]
    #修正优化模型
    responseStrategy = optimizationModel.responseModel(optimalDispatch, case, peak=peak, amount=amount,mode='E')
    #求解
    print('电需求响应' + str(100*e) + '%：'+ str(value(responseStrategy.objective)))
    #获取结果
    res_result = optimizationModel.retriveResult(microgrid_data,case,responseStrategy)
    res_result.to_excel(writer,sheet_name = str(e))
writer.save()
'''Heat DR'''
writer = pd.ExcelWriter('DayAhead Heat DR.xlsx')
peak = range(72, 76)
(max_model,max_amount) = optimizationModel.getMaxAmount(optimalDispatch,case,peak=peak,amount = [3000]*len(peak),mode='H')
for e in [1,0.8,0.6,0.4,0.2]:
    amount = [e * max_amount[t] for t in range(max_amount.__len__())]
    # 修正优化模型并求解
    responseStrategy = optimizationModel.responseModel(optimalDispatch, case, peak=peak, amount=amount,mode='H')
    print('热需求响应' + str(100 * e) + '%：' + str(value(responseStrategy.objective)))
    # 获取结果
    res_result = optimizationModel.retriveResult(microgrid_data, case, responseStrategy)
    res_result.to_excel(writer, sheet_name=str(e))
writer.save()