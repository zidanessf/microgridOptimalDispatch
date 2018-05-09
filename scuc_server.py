from importFile import *
import pandas as pd
import optimizationModel,microgridStructure
import scuc_pb2,scuc_pb2_grpc,grpc,time
from concurrent import futures
class SCUC(scuc_pb2_grpc.SCUCServicer):
    def GetSCUCResults(self,request,context):
        ppc = request.ppcfile
        tempfile = open('temp_ppc.py',mode='wb')
        tempfile.write(ppc)
        tempfile.close()
        case = getattr(importlib.import_module("temp_ppc"), request.ppc_method_name)
        ppc = case()
        graph_case = microgridStructure.ppc2graph(ppc)
        case_PS = microgridStructure.MicrogridCase_Graph(graph=graph_case, NumOfTime=96)
        load = pd.DataFrame()
        for load_one in request.load_all:
            bus = str(load_one.bus)
            load[bus] = pd.Series([x for x in load_one.load_of_a_day])
        mdl = optimizationModel.DayAheadModel(load,case_PS)
        solver = SolverManagerFactory('neos')
        res = solver.solve(mdl,opt=SolverFactory('cplex'))
        status = [scuc_pb2.STATUS(gen = int(gen.split('_')[1]),
                                     status_of_a_day = [int(value(mdl.gt_state[gen,t])) for t in mdl.T])\
                     for gen in mdl.N_gt]
        power = [scuc_pb2.POWER(gen = int(gen.split('_')[1]),
                                power_of_a_day = [value(mdl.gt_power[gen,t]) for t in mdl.T])\
                     for gen in mdl.N_gt]
        case_PS.update(mdl)
        case_PS.DCPowerFlow()
        power_flow = list()
        for branch in case_PS.graph.edges():
            nf = branch[0]
            nt = branch[1]
            power_flow.append(scuc_pb2.POWER_FLOW(
                f_bus = nf,t_bus = nt,power_flow_of_a_day=case_PS.graph.edge[nf][nt]['Power_Flow']))

        return scuc_pb2.SCUCoutput(termination_condition=str(res.solver.termination_condition),
                                   optimal_value=value(mdl.objective),
                                   status=status, power=power)
_ONE_DAY_IN_SECONDS = 60 * 60 * 24
def serve():
  server = grpc.server(futures.ThreadPoolExecutor(max_workers=10))
  scuc_pb2_grpc.add_SCUCServicer_to_server(SCUC(), server)
  server.add_insecure_port('[::]:44099')
  server.start()
  print('server started in port 44099')
  logging.info('server started in port 44099')

  try:
    while True:
      time.sleep(_ONE_DAY_IN_SECONDS)
  except KeyboardInterrupt:
    server.stop(0)
if __name__ == '__main__':
    serve()

# '''Initialize a special case of microgrid'''
# case = microgridStructure.case_PS
# '''Load input data'''
# microgrid_data = pd.read_excel('input_PS.xlsx')
# '''Construct base model'''
# optimalDispatch = optimizationModel.DayAheadModel(microgrid_data,case,range(96))
# solver = SolverFactory('cplex')
# res = solver.solve(optimalDispatch,tee=True)
