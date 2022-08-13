import math

from ac_metro.utility.abstractUtility import AbstractUtility

class BoundingBoxUtility(AbstractUtility):
    '''
    Create a bounding box of all stations in the input (with possible padding on ALL sides)
    '''

    @staticmethod
    def execute(map, options=None):

        if "padding" in options:
            pad = options["padding"]
        else:
            pad = 0

        min_x, min_y, max_x, max_y = math.inf, math.inf, -math.inf, -math.inf
        for i, (node, data) in enumerate(map.graph.nodes(data=True)):
            station = data["station"]
            min_x = min(min_x, station.x)
            max_x = max(max_x, station.x)
            min_y = min(min_y, station.y)
            max_y = max(max_y, station.y)

        return min_x - pad, min_y - pad, (max_x - min_x) + 2 * pad, (max_y - min_y) + 2 * pad
