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
optimalDispatch = optimizationModel.ConstructModel(microgrid_data,case)
'''Solve the base model'''
xfrm = TransformationFactory('gdp.chull')
xfrm.apply_to(optimalDispatch)
solver = SolverFactory('glpk')
solver.solve(optimalDispatch)
'''Retrieve the result'''
result = optimizationModel.retriveResult(microgrid_data,case,optimalDispatch)
result.to_excel('output.xlsx')
plt.plot(result)
plt.show()

'''Interaction process'''
'''Electric DR'''
'''
writer = pd.ExcelWriter('Electric DR.xlsx')
for e in [0.2,0.4,0.6,0.8,1]:
    amount = e * 3000
    responseStrategy = optimizationModel.responseModel(optimalDispatch, case, peak=range(72, 76), amount=amount,mode='E')
    xfrm = TransformationFactory('gdp.chull')
    xfrm.apply_to(responseStrategy)
    solver = SolverFactory('glpk')
    solver.solve(responseStrategy)
    res_result = optimizationModel.retriveResult(microgrid_data,case,responseStrategy)
    res_result.to_excel(writer,sheet_name = str(e))
writer.save()
'''
'''Heat DR'''
writer = pd.ExcelWriter('Heat DR.xlsx')
for e in [0.2,0.4,0.6,0.8,1]:
    amount = e * 3000
    responseStrategy = optimizationModel.responseModel(optimalDispatch, case, peak=range(72, 76), amount=amount,mode='H')
    xfrm = TransformationFactory('gdp.chull')
    xfrm.apply_to(responseStrategy)
    solver = SolverFactory('glpk')
    solver.solve(responseStrategy)
    res_result = optimizationModel.retriveResult(microgrid_data, case, responseStrategy)
    res_result.to_excel(writer, sheet_name=str(e))
writer.save()