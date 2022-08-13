import math

import numpy as np

from ac_metro.utility.abstractUtility import AbstractUtility


class CircularVertexOrderUtility(AbstractUtility):
    '''
    Compute a circular order at every station (returned as a dict), based on the current station positions (NOT taking controlpoints into account)
    '''

    @staticmethod
    def execute(map, options=None):
        circ_order = {}
        for id1, data in map.graph.nodes(data=True):

            node1 = data["station"]
            neighbours = []

            for id2 in map.graph.neighbors(id1):

                node2 = map.graph.nodes[id2]["station"]
                dX = node2.x - node1.x
                dY = node2.y - node1.y

                angle = np.arctan2(-dY, dX) * 180 / np.pi  # Important invert y axis

                if angle < 0:
                    angle = 360 + angle

                neighbours.append((id2, angle))

            if len(neighbours) < options["min_degree"]:
                continue

            neighbours.sort(key=lambda tup: tup[1])
            circ_order[id1] = [i[0] for i in neighbours]

        return circ_order