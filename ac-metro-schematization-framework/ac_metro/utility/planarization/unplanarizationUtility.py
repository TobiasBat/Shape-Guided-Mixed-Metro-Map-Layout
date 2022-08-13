from ac_metro.map.map import Map
from ac_metro.utility.abstractUtility import AbstractUtility
from ac_metro.utility.planarization.planarizationUtility import PlanarizationUtility


class UnPlanarizationUtility(AbstractUtility):
    '''
    Convert all added nodes back to intersections
    '''

    UTIL_NAME = 'UnPlanarizationUtility'

    @staticmethod
    def execute(map: Map, options=None):
        nodes = list(map.graph.nodes(data=True))

        for (node, nodeData) in nodes:
            if 'station' in nodeData:
                if nodeData['station'].hasContextualInformation(PlanarizationUtility.UTIL_NAME()):
                    map.removeStationKeepLines(node)

        return

