from ac_metro.utility.abstractUtility import AbstractUtility
from math import sin, cos

class MapTranslater(AbstractUtility):
    '''
    Use a dictionary with identifier-color pairs to color a metro map.
    '''

    @staticmethod
    def execute(map, options=None):
        if options == None:
            return

        if "translate_vector" not in options:
            print("No vector provided to TranslateUtility")
            return
        else:
            vector = options["translate_vector"]
            v_x = vector[0]
            v_y = vector[1]

        for i, (node, data) in enumerate(map.graph.nodes(data=True)):
            station = data["station"]
            # station.x, station.y = station.x + v_x, station.y + v_y
            station.x += v_x
            station.y += v_y
