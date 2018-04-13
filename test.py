from importFile import *
from pyomo.bilevel import *
from pyomo.mpec import *
m = ConcreteModel()
m.x = Var(bounds=(-1,1))
m.y = Var(bounds=(-1,1))
m.o = Objective(expr=m.x * m.x + m.y * m.y)
m.c = Constraint(expr=m.x * m.y == 1)
m.c1 = Constraint(expr=m.x == 1)
m.c2 = Constraint(expr=m.y == 1)
# solver = SolverFactory('gurobi')
# res = solver.solve(m)
# print(res)