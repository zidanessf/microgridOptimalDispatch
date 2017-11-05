from importFile import *
from microgrid_Model import *
import pandas as pd
import matplotlib.pyplot as plt
import optimizationModel,microgridStructure
'''Initialize a special case of microgrid'''
case = microgridStructure.MicrogridCase()
'''Load input data'''
microgrid_data = pd.read_excel('input.xlsx')
'''Construct base model'''
optimalDispatch = optimizationModel.DayAheadModel(microgrid_data,case)

'''Solve the base model'''

xfrm = TransformationFactory('gdp.chull')
xfrm.apply_to(optimalDispatch)
solver = SolverFactory('glpk')
solver.solve(optimalDispatch)
'''Retrieve the result'''
result = optimizationModel.retriveResult(microgrid_data,case,optimalDispatch)
result.to_excel('output.xlsx')
result.plot()
'''Interaction process'''
'''Electric DR'''
writer = pd.ExcelWriter('DayAhead Electric DR.xlsx')
peak = range(72, 76)
(max_model,max_amount) = optimizationModel.getMaxAmount(optimalDispatch,case,peak=peak,amount = [3000]*len(peak),mode='E')
max_result = optimizationModel.retriveResult(microgrid_data,case,max_model)
max_result.to_excel(writer,sheet_name = '1')
for e in [0.8,0.6,0.4,0.2]:
    amount = [e * max_amount[t] for t in range(max_amount.__len__())]
    #修正优化模型
    responseStrategy = optimizationModel.responseModel(optimalDispatch, case, peak=peak, amount=amount,mode='E')
    #求解
    xfrm = TransformationFactory('gdp.chull')
    xfrm.apply_to(responseStrategy)
    solver = SolverFactory('glpk')
    solver.solve(responseStrategy)
    #获取结果
    res_result = optimizationModel.retriveResult(microgrid_data,case,responseStrategy)
    res_result.to_excel(writer,sheet_name = str(e))
writer.save()
'''Heat DR'''
writer = pd.ExcelWriter('DayAhead Heat DR.xlsx')
peak = range(72, 76)
(max_model,max_amount) = optimizationModel.getMaxAmount(optimalDispatch,case,peak=peak,amount = [3000]*len(peak),mode='H')
max_result = optimizationModel.retriveResult(microgrid_data,case,max_model)
max_result.to_excel(writer,sheet_name = '1')
for e in [0.8,0.6,0.4,0.2]:
    amount = [e * max_amount[t] for t in range(max_amount.__len__())]
    # 修正优化模型并求解
    responseStrategy = optimizationModel.responseModel(optimalDispatch, case, peak=peak, amount=amount,mode='H')
    # 获取结果
    res_result = optimizationModel.retriveResult(microgrid_data, case, responseStrategy)
    res_result.to_excel(writer, sheet_name=str(e))
writer.save()