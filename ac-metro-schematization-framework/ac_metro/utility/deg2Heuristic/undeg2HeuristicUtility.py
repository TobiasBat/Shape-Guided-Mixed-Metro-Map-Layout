import math

from ac_metro.map.map import Map, Station
from ac_metro.utility.abstractUtility import AbstractUtility
from ac_metro.utility.deg2Heuristic.deg2HeuristicUtility import Deg2HeuristicUtility


class UnDeg2HeuristicUtility(AbstractUtility):
    '''
    Convert all added edges back to the original edges
    '''

    UTIL_NAME = 'UnDeg2HeuristicUtility'

    @staticmethod
    def execute(Map: Map, options=None):

        nodes = list(Map.graph.nodes(data=True))

        preserve_c = "preserve_control_points" in options and options["preserve_control_points"]

        for (node, nodeData) in nodes:
            if 'station' in nodeData:
                if nodeData['station'].hasContextualInformation(Deg2HeuristicUtility.UTIL_NAME()):
                    Map.convertStationToControlPoint(node)

        edges = list(Map.getEdges(None, True))

        for (startId, endId, edgeData) in edges:

            if startId != edgeData['edge'].source:
                temp = startId
                startId = endId
                endId = temp

            if 'edge' in edgeData:

                if preserve_c:
                    s_st = Map.getStation(edgeData['edge'].source)
                    t_st = Map.getStation(edgeData['edge'].target)
                    edge_geometry = [(s_st.x, s_st.y)] + edgeData['edge'].controlPoints + [(t_st.x, t_st.y)]
                    edge_proportions = [0]
                    edge_length = 0
                    for i in range(1, len(edge_geometry)):
                        segment_length = math.dist(edge_geometry[i - 1], edge_geometry[i])
                        edge_length += segment_length
                        edge_proportions.append(edge_length)

                if edgeData['edge'].hasContextualInformation(Deg2HeuristicUtility.UTIL_NAME()):
                    edgeInformation = edgeData['edge'].getContextualInformation(Deg2HeuristicUtility.UTIL_NAME())
                    if edgeInformation.numberOfReplacedEdges > 1 and edgeInformation.removedEdges is not None:
                        lines = list(edgeData['edge'].getLines())
                        for line in lines:
                            Map.removeConnection(startId, endId, line)

                        startData = Map.graph.nodes[startId]['station']
                        endData = Map.graph.nodes[endId]['station']
                        distanceX = endData.x - startData.x
                        distanceY = endData.y - startData.y

                        removedEdges = edgeInformation.removedEdges
                        if startId != removedEdges[0]['edge'].source and startId != removedEdges[0]['edge'].target:
                            removedEdges = list(reversed(removedEdges))

                        amountOfEdges = len(removedEdges)

                        previousStationData = startData
                        count = 0

                        if preserve_c:
                            c_point_step = edge_length / amountOfEdges
                            current_segment = 0

                        for removedEdge in removedEdges:
                            newStationId = removedEdge['edge'].target
                            if previousStationData.id == newStationId:
                                newStationId = removedEdge['edge'].source

                            if newStationId == endData.id:
                                if preserve_c:
                                    c_points = []
                                    current_segment_temp = current_segment
                                    while current_segment_temp < len(edge_geometry) - 1:
                                        c_points.append(edge_geometry[current_segment_temp])
                                        current_segment_temp += 1
                                for line in removedEdge['edge'].getLines():
                                    Map.addConnection(previousStationData.id, newStationId, line, controlPoints=c_points)
                                break

                            count += 1

                            # After this we know that the station will end up between the current_segment-1 and current_segment
                            if preserve_c:
                                c_points = []
                                while count * c_point_step >= edge_proportions[current_segment]:
                                    if current_segment > 0:
                                        c_points.append(edge_geometry[current_segment])
                                    current_segment += 1
                                length_on_current_segment = count * c_point_step - edge_proportions[current_segment - 1]
                                length_current_segment = math.dist(edge_geometry[current_segment - 1],
                                                                   edge_geometry[current_segment])
                                alpha = length_on_current_segment / length_current_segment
                                x_pos = (edge_geometry[current_segment - 1][0] * (1 - alpha)) + (
                                            edge_geometry[current_segment][0] * alpha)
                                y_pos = (edge_geometry[current_segment - 1][1] * (1 - alpha)) + (
                                            edge_geometry[current_segment][1] * alpha)
                            else:
                                x_pos = startData.x + (distanceX / (amountOfEdges)) * count
                                y_pos = startData.y + (distanceY / (amountOfEdges)) * count

                            newStationData = Station(
                                newStationId,
                                newStationId,
                                x_pos,
                                y_pos)
                            Map.addStation(newStationData)

                            for line in removedEdge['edge'].getLines():
                                Map.addConnection(previousStationData.id, newStationData.id, line,
                                                  controlPoints=c_points)

                            previousStationData = newStationData

        return
