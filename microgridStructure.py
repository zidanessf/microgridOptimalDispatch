from microgrid_Model import *
import networkx as nx
'''Initialize a special case of microgrid'''
#from neo4j.v1 import GraphDatabase

#driver = GraphDatabase.driver("bolt://localhost:7687", auth=("zidanessf", "123456"))
#graph = nx.Graph()
def get_node(tx, name, graph):
    return tx.run("MATCH (n:$name) RETURN n",name = name)


def get_branch(tx, name):
    for record in tx.run("MATCH (a:Person)-[:KNOWS]->(friend) WHERE a.name = $name "
                         "RETURN friend.name ORDER BY friend.name", name=name):
        print(record["friend.name"])

class MicrogridCase:
    def __init__(self):
        microgrid_device = dict()
        microgrid_device['PV_1'] = PV()
        #microgrid_device['ES_1'] = electricStorage()
        #microgrid_device['ABSC_1'] = absorptionChiller()
        microgrid_device['Boiler_1'] = boiler(Pmax=100)
        #microgrid_device['CS_1'] = coldStorage()
        microgrid_device['AC_1'] = airConditioner(Pmax=80,COP=0.9)
        microgrid_device['GT_1'] = gasTurbine(Pmax=500,efficiency = 0.33)
        microgrid_device['GT_2'] = gasTurbine(Pmax=500, efficiency=0.5)
        #microgrid_device['DR_Heat_Load'] = DRHeatLoad()
        microgrid_device['ut'] = utility()
        microgrid_device['inv'] = inverter(maxP=1000)
        microgrid_device['HS_1'] = heatStorage(capacity=1350,Hmax_in=150,Hmax_out=150)
        self.device = microgrid_device
        self.NumOfTime = 96
    def getKey(self,type):
        return [key for key in self.device.keys() if isinstance(self.device[key], type)]
