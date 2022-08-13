import math

import networkx as nx

from ac_metro.importing.abstractImport import AbstractImport


class ShapeImport(AbstractImport):
    '''
    Load a shape graph from an Graphml file.
    '''

    @staticmethod
    def load(options: dict):

        print(f"Loading Shape with cost {options['shape_cost']}")

        if options['filepath'] is None:
            raise FileNotFoundError

        g_in = nx.read_graphml(options['filepath'])

        if 'x_attr' in options:
            x_attr = options['x_attr']
        else:
            x_attr = 'x'

        if 'y_attr' in options:
            y_attr = options['y_attr']
        else:
            y_attr = 'y'

        if 'id_start' in options:
            id_counter = options['id_start']
        else:
            id_counter = 0

        shape_graph = nx.Graph()
        id_map = {}

        for v, data in g_in.nodes(data=True):
            x = float(data[x_attr])
            y = float(data[y_attr])
            id_map[v] = id_counter
            shape_graph.add_node(id_counter, x=x, y=y)
            id_counter += 1

        for u, v, data in g_in.edges(data=True):
            shape_graph.add_edge(id_map[u], id_map[v], crossing=[],
                                 weight=options["shape_cost"] if "shape_cost" in options else 1)

        return shape_graph
