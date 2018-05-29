from microgrid_Model import *
import networkx as nx
from scipy.sparse import dok_matrix
from scipy.linalg import pinv
from importFile import *
'''Initialize a special case of microgrid'''
class MicrogridCase_Simple:
    def __init__(self,device,NumOfTime,graph=None):
        self.device = device
        self.type = 'Simple'
        self.NumOfTime = NumOfTime
    def getKey(self,type):
        return [key for key in self.device.keys() if isinstance(self.device[key], type)]
    def SOCUpdate(self, plan, nowtime):
        N_es = self.getKey(electricStorage)
        for es in N_es:
            self.device[es].SOCnow = plan[es + '电池电量'].loc[nowtime] / self.device[es].capacity
class MicrogridCase_Graph:
    def __init__(self,graph,NumOfTime):
        ns = graph.nodes()
        # for i in range(len(ns)):
        #     graph.node[ns[i]].update({'index':i})
        self.graph = graph
        device = dict()
        for node in graph.nodes():
            device.update(graph.node[node]['device'])
        self.device = device
        self.type = 'Graph'
        self.NumOfTime = NumOfTime
        self.Adj = nx.adjacency_matrix(graph)
        ns = graph.nodes()
        B = dok_matrix((len(ns),len(ns)))
        for l in graph.edges():
            B[l[0],l[1]] = - 1/graph.edge[l[0]][l[1]]['X']
        B = B + B.transpose()
        for i in range(len(ns)):
            B[i,i] = - sum(B[i,j] for j in range(len(ns)))
        self.B = B
        self.B_INV = pinv(B.todense())
    def getKey(self, type):
        return [key for key in self.device.keys() if isinstance(self.device[key], type)]
    def DCPowerFlow(self):
        for branch in self.graph.edges():
            nf = branch[0]
            nt = branch[1]
            X = self.graph.edge[nf][nt]['X']
            PF = list()
            for t in self.T:
                PF.append(1 / X * (
                sum(self.B_INV.item(nf, nb) * self.graph.node[nb]['P_inj'][t] for nb in self.graph.nodes()) - sum(
                    self.B_INV.item(nt, nb) * self.graph.node[nb]['P_inj'][t] for nb in self.graph.nodes())))
            self.graph.edge[nf][nt].update({'Power_Flow' : PF})
    def update(self,mdl):
        def Power_Injection(mdl,graph, nb, t):
            Temp = 0.0
            for key, dev in graph.node[nb]['device'].items():
                if isinstance(dev, gasTurbine):
                    Temp += mdl.sub.gt_power[key, t]
                    dev.result.append(value(mdl.sub.gt_power[key, t]))
                if isinstance(dev, PV):
                    Temp += mdl.wp[key, t]
                    dev.result.append(value(mdl.wp[key, t]))
            return value(Temp)
        for bus in self.graph.nodes():
            self.graph.node[bus].update({'P_inj' : [Power_Injection(mdl,self.graph,bus,t) for t in mdl.sub.T]})
        self.T = mdl.sub.T

device_IES = {
    'PV_1' : PV(),
    'ES_1' : electricStorage(),
    'ABSC_1' : absorptionChiller(),
    'BOL_1' : boiler(),
    'CS_1' : coldStorage(),
    'AC_1' : airConditioner(),
    'GT_1' : gasTurbine(),
    'DR_Heat_Load' : DRHeatLoad(),
    'ut' : utility(),
    'inv' : inverter()
}
case_IES = MicrogridCase_Simple(device=device_IES, NumOfTime=96)
# graph_PS = nx.Graph()
# graph_PS.add_nodes_from([0,1,2,3,4])
# graph_PS.add_edges_from([(0,1),(0,4),(1,2),(2,3),(3,4),(0,3)])
# graph_PS.node[0].update({
#     'ID' : 'A',
#     'device' : {
#         'Park City' : gasTurbine(Pmax=170,Pmin=10,Cost=15),
#         'Alta' : PV()
#     }
# })
# graph_PS.node[1].update({
#     'ID' : 'B',
#     'device': {}
#     })
# graph_PS.node[2].update({
#     'ID' : 'C',
#     'device':{
#         'Solitude' : gasTurbine(Pmax=520,Pmin=10,Cost=30)
#     }
# })
# graph_PS.node[3].update({
#     'ID' : 'D',
#     'device':{
#         'Sundance' : gasTurbine(Pmax=200,Pmin=10,Cost=40)
#     }
# })
# graph_PS.node[4].update({
#     'ID' : 'E',
#     'device':{
#         'Brighton':gasTurbine(Pmax=600,Pmin=10,Cost=20)
#     }
# })
# graph_PS.edge[0][1].update({
#     'R' : 0.281,
#     'X' : 2.81,
#     'Limit' : 400
# })
# graph_PS.edge[0][3].update({
#     'R' : 0.304,
#     'X' : 3.04,
#     'Limit' : None
# })
# graph_PS.edge[0][4].update({
#     'R' : 0.064,
#     'X' : 0.64,
#     'Limit' : None
# })
# graph_PS.edge[1][2].update({
#     'R' : 0.108,
#     'X' : 1.08,
#     'Limit' : None
# })
# graph_PS.edge[2][3].update({
#     'R' : 0.297,
#     'X' : 2.97,
#     'Limit' : None
# })
# graph_PS.edge[3][4].update({
#     'R' : 0.297,
#     'X' : 2.97,
#     'Limit' : 240
# })
# 赶工啦，临时写的一个PYPOWER ---- GRAPH 转换器
from pypowercase.case39 import case39
import numpy as np
np.random.seed(96)
def ppc2graph(ppc):
    graph_PS = nx.Graph()
    ppc_nodes = ppc['bus'][:,0].astype(np.int) - 1
    graph_PS.add_nodes_from(ppc_nodes)
    for bus in graph_PS.nodes():
        graph_PS.node[bus].update({
            'device':{},
            'Load' :  (0.3*np.random.rand(96) + 0.85)*[ppc['bus'][bus,2]]
        })
    gen = ppc['gen']
    gencost = ppc['gencost']
    for row in gen:
        bus = row[0].astype(np.int) - 1
        graph_PS.node[bus]['device'].update({'GT_'+str(bus):gasTurbine(Pmax=row[8],Pmin=row[9],
                                                                       Cost=0.5*bus)})
    for row in ppc['branch']:
        nf = row[0].astype(np.int) - 1
        nt = row[1].astype(np.int) - 1
        graph_PS.add_edge(nf,nt,{'R':row[2],'X':row[3],'Limit':row[5]})

    return graph_PS

ppc = case39()
graph_case39 = ppc2graph(ppc)
graph_case39.node[1]['device'].update({'WT':PV()})
case_PS = MicrogridCase_Graph(graph=graph_case39,NumOfTime=96)
