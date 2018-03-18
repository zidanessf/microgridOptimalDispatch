from importFile import *
from microgrid_Model import *
import pandas as pd
import matplotlib.pyplot as plt
import optimizationModel,microgridStructure
'''Initialize a special case of microgrid'''
case = microgridStructure.MicrogridCase_Simple()
'''Load input data'''
microgrid_data = pd.read_excel('input.xlsx')
Econtract = pd.read_excel('DayAhead Electric DR.xlsx',sheetname='1')
Hcontract = []
case.SOCUpdate(plan = Econtract, nowtime = 4)
input_data = {
    "microgrid_data" : microgrid_data,
    "case" : case,
    "refE" : Econtract['购电功率'],
    "refH" : [],
    "refSS" : {
        'GT_1':[1 if p > 0.1 else 0 for p in Econtract['GT_1机组出力']],
        'Boiler_1':[1 if p > 0.1 else 0 for p in Econtract['Boiler_1燃气锅炉热功率']]
    },
    "peak" : range(72,76),
    "T_range" : range(4,96)
}
'''Construct base model'''
optimalDispatch = optimizationModel.DayInModel(**input_data)
'''Solve the base model'''
xfrm = TransformationFactory('gdp.chull')
xfrm.apply_to(optimalDispatch)
solver = SolverFactory('glpk')
solver.solve(optimalDispatch)
# optimalDispatch.P_0 = [value(optimalDispatch.utility_power[t]) for t in optimalDispatch.T]
# optimalDispatch.H_0 = [value(optimalDispatch.buy_heat[t]) for t in optimalDispatch.T]
print('自趋优：'+str(value(optimalDispatch.objective)))
'''Retrieve the result'''
result = optimizationModel.retriveResult(microgrid_data,case,optimalDispatch)
result.to_excel('dayin.xlsx')
optimizationModel.extendedResult(result)
'''Interaction process'''
'''Electric DR'''
writer = pd.ExcelWriter('DayIn Electric DR.xlsx')
peak = range(82, 86)
(max_model,max_amount) = optimizationModel.getMaxAmount(optimalDispatch,case,peak=peak,amount = [3000]*len(peak),mode='E')
for e in [1,0.8,0.6,0.4,0.2]:
    amount = [e * max_amount[t] for t in range(max_amount.__len__())]
    #修正优化模型
    responseStrategy = optimizationModel.responseModel(optimalDispatch, case, peak=peak, amount = amount,mode='E')
    #求解
    xfrm = TransformationFactory('gdp.chull')
    xfrm.apply_to(responseStrategy)
    solver = SolverFactory('glpk')
    solver.solve(responseStrategy)
    print('电需求响应' + str(100 * e) + '%：' + str(value(responseStrategy.objective)))
    #获取结果
    res_result = optimizationModel.retriveResult(microgrid_data,case,responseStrategy)
    '''电关口功率曲线画图，并保存excel'''
    plt.plot(res_result['购电功率'],label =  (str(100 * e) + '%电需求响应',))
    plt.legend()
    res_result.to_excel(writer,sheet_name = str(e))
writer.save()
plt.show()
'''Heat DR'''
writer = pd.ExcelWriter('DayIn Heat DR.xlsx')
peak = range(82, 86)
(max_model,max_amount) = optimizationModel.getMaxAmount(optimalDispatch,case,peak=peak,amount = [3000]*len(peak),mode='H')
plt.figure(2)
for e in [1,0.8,0.6,0.4,0.2]:
    amount = [e * max_amount[t] for t in range(max_amount.__len__())]
    # 修正优化模型并求解
    responseStrategy = optimizationModel.responseModel(optimalDispatch, case, peak=peak, amount=amount,mode='H')
    print('热需求响应' + str(100 * e) + '%：' + str(value(responseStrategy.objective)))
    # 获取结果
    res_result = optimizationModel.retriveResult(microgrid_data, case, responseStrategy)
    '''热关口功率曲线画图，并保存excel'''
    plt.plot(res_result['购热功率'],label=str(100 * e) + '%热需求响应')
    plt.legend()
    res_result.to_excel(writer, sheet_name=str(e))
writer.save()
plt.show()