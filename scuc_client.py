import grpc,scuc_pb2_grpc,scuc_pb2,case39,microgridStructure
ppc = open('./cases/case39.py')
channel = grpc.insecure_channel('localhost:44099')
stub = scuc_pb2_grpc.SCUCStub(channel)
ppc_bytes = bytes(ppc.read(),encoding='utf-8')
ppc.close()
ppc_graph = microgridStructure.ppc2graph(case39.case39())
load_all = [scuc_pb2.LOAD(bus=bus,load_of_a_day=ppc_graph.node[bus]['Load']) for bus in ppc_graph.nodes()]
response = stub.GetSCUCResults(scuc_pb2.SCUCinput(ppc_method_name = 'case39',ppcfile=ppc_bytes,load_all=load_all))
print(response.termination_condition)
