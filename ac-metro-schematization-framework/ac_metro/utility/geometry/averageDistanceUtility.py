import math

from ac_metro.utility.abstractUtility import AbstractUtility


class AverageDistanceUtility(AbstractUtility):
    '''
    Compute the average distance between two stations
    '''

    @staticmethod
    def execute(map, options=None):
        d_total = 0
        data = map.graph.nodes(data=True)
        for u, v in map.graph.edges():
            # print(f"for this edge ({u}, {v}) we get the x coordinates\nu_x: {map.graph.nodes(data=True)[u]['station'].x}\nu_y: {map.graph.nodes[u]['station'].y}\nv_x: {map.graph.nodes[v]['station'].x}\nv_y: {map.graph.nodes[v]['station'].y}")
            x_u, y_u = data[u]['station'].x, data[u]['station'].y
            x_v, y_v = data[v]['station'].x, data[v]['station'].y
            d_total += math.sqrt(math.pow(x_u - x_v, 2) + math.pow(y_u - y_v, 2))
        return d_total / map.graph.number_of_edges()
