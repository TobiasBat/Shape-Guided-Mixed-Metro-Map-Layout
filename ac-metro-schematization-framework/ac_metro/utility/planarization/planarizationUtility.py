from shapely.geometry import LineString

from ac_metro.map.map import Map, Edge
from ac_metro.utility.abstractUtility import AbstractUtility

class AnnotatedEdge:
    def __init__(self, s1, s2, edge: LineString):
        self.s1 = s1
        self.s2 = s2
        self.edge = edge

    def __str__(self):
        return f'({self.s1} <-> {self.s2})'

    def __repr__(self):
        return f'({self.s1} <-> {self.s2})'

def getEndpointPositions(edge: Edge, map: Map) -> AnnotatedEdge:
    (start, end, _) = edge
    startNode = map.graph.nodes[start]
    endNode = map.graph.nodes[end]
    return AnnotatedEdge(
        startNode['station'].id,
        endNode['station'].id,
        LineString([
            (startNode['station'].x, startNode['station'].y),
            (endNode['station'].x, endNode['station'].y)
        ])
    )

class PlanarizationUtility(AbstractUtility):
    '''
    Convert intersections to nodes
    '''

    @staticmethod
    def UTIL_NAME():
        return 'PlanarizationUtility'

    @staticmethod
    def execute(map: Map, options=None):
        foundIntersection = True
        count = 0
        while (foundIntersection):
            foundIntersection = False

            edges = map.getEdges(None, True)
            mappedEdges = []
            for edge in edges:
                mappedEdges.append(getEndpointPositions(edge, map))

            for i in range(len(mappedEdges)):
                for j in range(i + 1, len(mappedEdges)):
                    annotatedEdgeA = mappedEdges[i]
                    annotatedEdgeB = mappedEdges[j]
                    # print(f"Checking {annotatedEdgeA} and {annotatedEdgeB}")
                    if annotatedEdgeA.s1 == annotatedEdgeB.s1 \
                        or annotatedEdgeA.s1 == annotatedEdgeB.s2 \
                        or annotatedEdgeA.s2 == annotatedEdgeB.s1 \
                        or annotatedEdgeA.s2 == annotatedEdgeB.s2:
                        continue
                    # print(f"skipping {list(annotatedEdgeA.edge.coords)[0]} and {list(annotatedEdgeB.edge.coords)}?")
                    # if list(annotatedEdgeA.edge.coords)[0] == list(annotatedEdgeA.edge.coords)[1]\
                    #     or list(annotatedEdgeB.edge.coords)[0] == list(annotatedEdgeB.edge.coords)[1]:
                    #     continue

                    if annotatedEdgeA.edge.intersects(annotatedEdgeB.edge):
                        intersection = annotatedEdgeA.edge.intersection(annotatedEdgeB.edge)
                        if list(intersection.coords)[0] == list(annotatedEdgeA.edge.coords)[0] \
                            or list(intersection.coords)[0] == list(annotatedEdgeA.edge.coords)[1] \
                            or list(intersection.coords)[0] == list(annotatedEdgeB.edge.coords)[0] \
                            or list(intersection.coords)[0] == list(annotatedEdgeB.edge.coords)[1]:
                            # print(f"intersection at {list(intersection.coords)[0]}")
                            # print(f"coords: {list(annotatedEdgeA.edge.coords)[0]}")
                            # print(f"coords: {list(annotatedEdgeA.edge.coords)[1]}")
                            # print(f"coords: {list(annotatedEdgeB.edge.coords)[0]}")
                            # print(f"coords: {list(annotatedEdgeB.edge.coords)[1]}")
                            continue

                        print(f"Intersect: {annotatedEdgeA} and {annotatedEdgeB}")
                        map.convertIntersectionToStation(
                            annotatedEdgeA.s1,
                            annotatedEdgeA.s2,
                            annotatedEdgeB.s1,
                            annotatedEdgeB.s2,
                            f"intersection_{count}",
                            PlanarizationUtility.UTIL_NAME()
                        )
                        count += 1
                        foundIntersection = True
                        break
                if foundIntersection:
                    break


