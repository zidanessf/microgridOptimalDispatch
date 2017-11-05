from importFile import *
import copy
model = ConcreteModel()
model.x_1 = Var(within=NonNegativeReals)
model.x_2 = Var(within=NonNegativeReals)
model.obj = Objective(expr=5*(model.x_1 + 2*model.x_2))
model.con1 = Constraint(expr=3*model.x_1 + 4*model.x_2 >= 1)
model.con2 = Constraint(expr=2*model.x_1 + 5*model.x_2 >= 2)
'''
tmpmdl = copy.deepcopy(model)
tmpmdl.temp_obj = Objective(expr=model.x_1 + model.x_2)
model.obj.set_value(model.obj.expr + tmpmdl.temp_obj.expr)
del tmpmdl
'''

solver = SolverFactory('glpk')
solver.solve(model)