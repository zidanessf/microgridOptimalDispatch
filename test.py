from neo4j.v1 import GraphDatabase
import networkx as nx
driver = GraphDatabase.driver("bolt://localhost:7687", auth=("zidanessf", "123456"))
G = nx.Graph()
def get_node(tx,keyword,graph):
    mylist = list()
    for record in  tx.run("MATCH (n:{}) RETURN n",):
        graph.add_node(record._values[0].id,record._values[0].properties)
    return graph
s = driver.session()
result = s.read_transaction(get_node,G)
1