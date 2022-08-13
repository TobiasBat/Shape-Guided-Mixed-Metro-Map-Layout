import math

from ac_metro.map.map import Map, Station, Edge
from ac_metro.utility.abstractUtility import AbstractUtility

class NodeLink:
    def __init__(self, nodeId):
        self.nodeId = nodeId
        self.previous = None
        self.next = None

    def setPrevious(self, previous):
        self.previous = previous
        previous.next = self

    def setNext(self, next):
        self.next = next
        next.previous = self

    def __str__(self):
            return f'({ self.previous.nodeId if self.previous is not None else "" } <- {self.nodeId} -> { self.next.nodeId if self.next is not None else "" })'

    def __repr__(self):
        return f'({ self.previous.nodeId if self.previous is not None else "" } <- {self.nodeId} -> { self.next.nodeId if self.next is not None else "" })'


class Path:
    def __init__(self, startId, endId, innerStations: [Station], innerEdges: [Edge]):
        self._startId = startId
        self._endId = endId
        self._innerStations = innerStations
        self._innerEdges = innerEdges

    def getStartId(self):
        return self._startId

    def getEndId(self):
        return self._endId

    def getInnerStations(self):
        return self._innerStations

    def getInnerIDs(self):
        return [self._startId] + list(map((lambda station: station.id), self._innerStations)) + [self._endId]

    def getInnerEdges(self):
        return self._innerEdges

    def __str__(self):
        return f'({self._startId} ... {self._endId})'

    def __repr__(self):
        return f'({self._startId} ... {self._endId})'


class PathInformation:

    def __init__(self, path_array: [Path], path_edge_dictionary):
        self.path_array = path_array
        self.path_edge_dictionary = path_edge_dictionary

    def getPaths(self):
        return self.path_array

    def edgeIn(self, startId, endId):
        if (startId, endId) in self.path_edge_dictionary:
            return self.path_edge_dictionary[(startId, endId)]
        if (endId, startId) in self.path_edge_dictionary:
            return self.path_edge_dictionary[(endId, startId)]
        return None

def getPath(startId, endId, map: Map) -> Path:

    start = NodeLink(startId)
    end = NodeLink(endId)
    start.setNext(end)

    while True:
        neighborsOfStart = list(map.graph.neighbors(start.nodeId))
        if len(neighborsOfStart) == 2:
            newStartId = neighborsOfStart[0]
            if newStartId == start.next.nodeId:
                newStartId = neighborsOfStart[1]
            if newStartId == end.nodeId: # check for cycle
                break
            newStart = NodeLink(newStartId)
            start.setPrevious(newStart)
            start = newStart
        else:
            break

    while True:
        neighborsOfEnd = list(map.graph.neighbors(end.nodeId))
        if len(neighborsOfEnd) == 2:
            newEndId = neighborsOfEnd[0]
            if newEndId == end.previous.nodeId:
                newEndId = neighborsOfEnd[1]
            if newEndId == start.nodeId: # check for cycle
                break
            newEnd = NodeLink(newEndId)
            end.setNext(newEnd)
            end = newEnd
        else:
            break

    currentNodeLink = start
    innerStations = []
    innerEdges = []
    while currentNodeLink.next is not end:
        nextNodeLink = currentNodeLink.next
        innerStations.append(map.graph.nodes[nextNodeLink.nodeId]['station'])
        innerEdges.append(map.graph.get_edge_data(currentNodeLink.nodeId, nextNodeLink.nodeId))
        currentNodeLink = nextNodeLink
    innerEdges.append(map.graph.get_edge_data(currentNodeLink.nodeId, end.nodeId))

    return Path(start.nodeId, end.nodeId, innerStations, innerEdges)


def removeEdgeFromList(id1, id2, edgeList: [Edge]):
    for i in range(len(edgeList)):
        (startId, endId, _) = edgeList[i]
        if (startId == id1 and endId == id2) or (startId == id2 and endId == id1):
            edgeList.pop(i)
            return


class GetPathsUtility(AbstractUtility):
    '''
    Finds all paths denoted by inner vertices of degree 2, that start and end in the nearest vertices with degree other than two
    '''

    @staticmethod
    def UTIL_NAME():
        return 'GetPathsUtility'

    @staticmethod
    def execute(Map, options={}):

        path_array = []
        edge_in_dictionary = {}

        edges = list(Map.getEdges(None, True))
        count = 0
        while len(edges) > 0:
            (startId, endId, edgeData) = edges.pop()

            path = getPath(startId, endId, Map)

            if len(path.getInnerStations()) > 0:
                edge_in_dictionary[(path.getStartId(), path.getInnerStations()[0].id)] = count
                removeEdgeFromList(path.getStartId(), path.getInnerStations()[0].id, edges)

                edge_in_dictionary[(path.getInnerStations()[-1].id, path.getEndId())] = count
                removeEdgeFromList(path.getInnerStations()[-1].id, path.getEndId(), edges)
            else:
                edge_in_dictionary[(path.getStartId(), path.getEndId())] = count

            for i in range(1, len(path.getInnerStations())):
                edge_in_dictionary[(path.getInnerStations()[i - 1].id, path.getInnerStations()[i].id)] = count
                removeEdgeFromList(path.getInnerStations()[i - 1].id, path.getInnerStations()[i].id, edges)

            path_array.append(path)
            count += 1

        return PathInformation(path_array, edge_in_dictionary)


