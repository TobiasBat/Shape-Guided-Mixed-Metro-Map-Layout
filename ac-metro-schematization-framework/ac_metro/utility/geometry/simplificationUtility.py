from ac_metro.utility.abstractUtility import AbstractUtility
import simplification.cutil as simp
import networkx as nx


class Simplifier(AbstractUtility):

    @staticmethod
    def execute(map, options=None):
        graph = options["graph"] if "graph" in options else map.graph
        factor = options["factor"] if "factor" in options else 1
        id_start = options["id_start"] if "id_start" in options else 0

        coord_list = graph_to_coord_list(graph)
        simplified = simp.simplify_coords(coord_list, factor)

        simplified_graph = nx.Graph()
        ID = id_start
        for coordinate in simplified:
            simplified_graph.add_node(ID, x=coordinate[0], y=coordinate[1])
            if ID > id_start:
                simplified_graph.add_edge(ID - 1, ID, crossing=[])
            ID += 1

        simplified_graph.add_edge(ID - 1, id_start, crossing=[])
        return simplified_graph


# Returns a list of tuples of coordinates of all vertices in the grid
def graph_to_coord_list(graph):
    coordinates = []
    for _, data in graph.nodes(data=True):
        coordinates.append([data["x"], data["y"]])
    return coordinates
