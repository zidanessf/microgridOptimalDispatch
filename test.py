from importFile import *
from pyomo.bilevel import *
from pyomo.mpec import *
import requests
model = ConcreteModel()
model.x = Var(bounds=(1,2))
model.v = Var(bounds=(1,2))
model.t = Var(within=Binary)
model.c = Constraint(expr=model.x + model.v <= 3)
model.o = Objective(expr=2*model.x + 2*model.v - model.t)
solver = SolverFactory('glpk')
result = solver.solve(model)
print(value(model))

#总结：pyomo可以处理上层带有整数变量或互补条件，下层为线性模型的双层规划模型。通过bilevel.linear_mpec进行转化，并通过bilevel_blp_global求解器求解