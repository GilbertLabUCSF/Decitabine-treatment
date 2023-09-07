"""
https://github.com/abearab/FunGI/blob/main/search.py
"""

import pandas as pd
import igraph as ig
import networkx as nx


# """
#   - Creating a graph
# """
# my_graph = ig.Graph.TupleList(
#     df[['col1','col2']].itertuples(index=False), 
#     directed=True, weights=False
# )#, edge_attrs="weight")


def node_i(graph, node):
    return [i for i, x in enumerate(graph.vs['name']) if x == node][0]


def filter_graph_by_nodes_BFS(G, nodes, max_distance=1, verbose=True):
    """filter graph `G` by running BFS for given nodes
    """
    bfs = []
    for node in nodes:
        for distance in range(1, max_distance + 1):
            bfs = bfs + run_bfs(G, node, distance=distance)
    
    outG = G.vs.select(name_in=bfs + nodes).subgraph()
    
    if verbose: get_graph_features(outG)
    
    return outG


def run_bfs(graph, node, distance=1):
    """Run BFS – Breadth-First Traversal (or Search)
    This function is running BFS algorithm to create a subgraph for given 1) `igraph` object, 2) node id, and 3) distance from the node.
    The subgraph includes all nodes with certain distance (defualt 1) from the given node.
    """
    A = [edge.tuple for edge in graph.es]
    G = nx.Graph(A)
    bfs = nx.descendants_at_distance(G, source=node_i(graph, node), distance=distance)
    out = [graph.vs['name'][i] for i in list(bfs)]
    return out


def nodes2df(G):
    node_df = pd.DataFrame({attr: G.vs[attr] for attr in G.vertex_attributes()})
    return node_df


def get_graph_features(g):
    print(g.summary(), '\n')
    print("Number of vertices in the graph:", g.vcount())
    print("Number of edges in the graph", g.ecount())
    print("Is the graph directed:", g.is_directed())
    print("Maximum degree in the graph:", g.maxdegree())

    top_nodes = g.vs.select(_degree=g.maxdegree())["name"]
    print("Node name with Maximum degree:", top_nodes)


def plot_graph(G, geneset=None, geneset2=None, layout="kk", vs_label_size=6, vs_size=35, b1=500, b2=500):
    """visualising graph data
    """
    G.vs["color"] = ["lightgray" for vertex in G.vs] 

    if geneset: 
        for gene in geneset:
            G.vs[[i for i, x in enumerate(G.vs['name']) if x == gene][0]]['color'] = 'yellow'
            
    if geneset2: 
        for gene in geneset2:
            G.vs[[i for i, x in enumerate(G.vs['name']) if x == gene][0]]['color'] = "lightblue"

    return ig.plot(
        G,
        layout=G.layout(layout),
        vertex_label=G.vs["name"],
        vertex_color=G.vs["color"],
        vertex_label_size=vs_label_size, 
        vertex_size=vs_size,
        edge_arrow_size = vs_label_size / 12,
        edge_arrow_width = vs_label_size / 3,
        bbox=(b1, b2), margin=60,
    )