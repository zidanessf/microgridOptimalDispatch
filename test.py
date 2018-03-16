from importFile import *
from pyomo.bilevel import *
from pyomo.mpec import *
import requests
model = ConcreteModel()
model.x = Var(bounds=(1,2))
model.v = Var(bounds=(1,2))
model.t = Var(within=Binary)
model.sub = SubModel()
model.sub.y = Var(bounds=(-2,5))
model.sub.w = Var(bounds=(-1,5))
model.o = Objective(expr=2*model.x + 2*model.sub.y + 2*model.v - model.t)
model.c = Constraint(expr=model.x + model.v >= 0)
model.sub.o = Objective(expr=model.x + model.sub.w, sense=maximize)
model.sub.c = Constraint(expr=2.29 <= model.sub.y + model.sub.w <= 2.3)

def balance(mdl):
    return mdl.y <= mdl.w
model.sub.balance = Constraint(rule=balance)
xfrm = TransformationFactory('bilevel.linear_mpec')
xfrm.apply_to(model)
solver = SolverFactory('bilevel_blp_global')
result = solver.solve(model)
print(value(model.sub.y))

#总结：pyomo可以处理上层带有整数变量或互补条件，下层为线性模型的双层规划模型。通过bilevel.linear_mpec进行转化，并通过bilevel_blp_global求解器求解