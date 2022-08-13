import collections
import copy
import networkx as nx
import datetime
import math
from shapely.geometry import LineString, Point

class ContextualInformation:

    def __init__(self):
        self.information = {}

    def hasInformation(self, author):
        return author in self.information

    def getInformation(self, author):
        return self.information[author]

    def setInformation(self, author, information):
        self.information[author] = information

    def removeInformation(self, author):
        self.information.pop(author)

    def getAuthors(self):
        return self.information.keys()

class MapSnapshot:

    def __init__(self, map, name=''):
        self.map = map
        self.name = name
        self.time = datetime.datetime.now().isoformat()

class Map:

    def __init__(self):
        self.graph = nx.Graph()
        self.lines = dict()
        self.lineGraphs = dict()
        self.lineColors = collections.defaultdict(type("#000000"))
        self.snapshots = []

    def getLines(self):
        '''
        Convenient getter method for lines in the network. Primarily used in for loops.

        :return: (key, value) iterator over all lines in the map
        '''
        return self.lines.items()

    def getLineGraphs(self):
        '''
        Iterator over all lines in the network. Useful in cases where the lines are not strictly consistent of a linear order over the stations (e.g. branching).

        :return: (key, value) iterator over all lines and their induced subgraphs.
        '''
        return self.lineGraphs.items()

    def addStation(self, station):
        '''
        Function to add a station to the map. The station must be created before being added. Station can further be referenced by its ID in the graph.

        :station: station to be added to the map. Creation of the station is provided by the user.
        
        :return: None
        '''
        self.graph.add_node(station.id, station=station)

    def addStationsBetween(self, station1Id, station2Id, amountOfStations, ids, creatorName=None, stationInformation=None, edgeInformation=None):
        '''
        Inserts stations equidistantly along the given edge. This edge must not contain control points.

        :return: None
        '''
        if self.graph.has_edge(station1Id, station2Id):
            edge = self.graph.get_edge_data(station1Id, station2Id)['edge']

            edge_points = [(self.graph.nodes[station1Id]['station'].x, self.graph.nodes[station1Id]['station'].y)] \
                + edge.getControlPoints(station1Id) \
                + [(self.graph.nodes[station2Id]['station'].x, self.graph.nodes[station2Id]['station'].y)]
            line_string = LineString(edge_points)

            control_point_positions = list(map(lambda controlpoint: (line_string.project(Point(controlpoint), True), controlpoint), edge.getControlPoints(station1Id)))


            lines = list(edge.getLines())
            for line in lines:
                self.removeConnection(station1Id, station2Id, line)

            previousStationId = station1Id

            for i in range(amountOfStations):
                station = Station(
                    ids[i],
                    ids[i],
                    line_string.interpolate((i + 1) / (amountOfStations + 1), True).x,
                    line_string.interpolate((i + 1) / (amountOfStations + 1), True).y,
                    creatorName,
                    stationInformation[i] if isinstance(stationInformation, list) else None)
                self.addStation(station)

                for line in lines:
                    self.addConnection(previousStationId, station.id, line, creatorName, edgeInformation[i] if isinstance(edgeInformation, list) else None)
                control_points_to_add = []
                while len(control_point_positions) > 0:
                    next_control_point_position = control_point_positions[0]
                    if (i + 1) / (amountOfStations + 1) >= next_control_point_position[0]:
                        control_points_to_add.append(control_point_positions.pop(0)[1])
                    else:
                        break
                self.getEdge(previousStationId, station.id).setControlPoints(control_points_to_add, previousStationId)

                previousStationId = station.id

            for line in lines:
                self.addConnection(previousStationId, station2Id, line, creatorName, edgeInformation[-1] if isinstance(edgeInformation, list) else None)
            self.getEdge(previousStationId, station2Id).setControlPoints(list(map(lambda control_point: control_point[1], control_point_positions)), previousStationId)

    def removeStation(self, stationId):
        '''
        Function to remove a station from the map. The station must not be connected to any lines (or edges).

        :stationId: station to be removed.

        :return: None
        '''
        if self.graph.has_node(stationId):
            # any(True for _ in iterator)
            if len(list(self.graph.neighbors(stationId))) == 0: # can be made more efficient, only check whether iterator emits one value
                self.graph.remove_node(stationId)
            else:
                raise Exception(f'Station "{ stationId }" to be removed although it still has lines connected.')

    def removeStationKeepLines(self, stationId):
        '''
        Function to remove a station from the map and keep the lines going through. There must not be more than two neighbors that are connected by the same line.

        :stationId: station to be removed.

        :return: None
        '''
        if self.graph.has_node(stationId):
            lineConnections = {}
            controlPointsTo = {}
            neighbors = list(self.graph.neighbors(stationId))
            for neighbor in neighbors:
                edge = self.graph.get_edge_data(stationId, neighbor)['edge']
                controlPointsTo[neighbor] = edge.getControlPoints(stationId)
                lines = edge.getLines()
                for line in list(lines):
                    if line in lineConnections:
                        lineConnections[line].append(neighbor)
                    else:
                        lineConnections[line] = [neighbor]
                    self.removeConnection(stationId, neighbor, line)

            self.removeStation(stationId)

            for (line, stations) in lineConnections.items():
                if len(stations) == 2:
                    self.addConnection(stations[0], stations[1], line)
                    # performance can be improved by setting this information only when edge is first created
                    self.graph.get_edge_data(stations[0], stations[1])['edge'].setControlPoints(
                        list(reversed(controlPointsTo[stations[0]])) + controlPointsTo[stations[1]],
                        stations[0]
                    )
                elif len(stations) == 1:
                    if len(neighbors) == 2:
                        otherNeighbor = neighbors[0]
                        if otherNeighbor == stations[0]:
                            otherNeighbor = neighbors[1]
                        self.addConnection(stations[0], otherNeighbor, line)
                        self.graph.get_edge_data(stations[0], otherNeighbor)['edge'].setControlPoints(
                            list(reversed(controlPointsTo[stations[0]])) + controlPointsTo[otherNeighbor],
                            stations[0]
                        )
                    else:
                        raise Exception(f'Only 1 out of more than 2 neighbors connected by the line "{ line }".')
                else:
                    raise Exception(f'More than 2 neighbors connected by the line "{ line }".')

        return

    def convertStationToControlPoint(self, stationId):
        '''
        Function to remove a station of degree two from the map and replaces it with a control point on the lines going through.

        :stationId: station to be converted.

        :return: None
        '''
        if self.graph.has_node(stationId):
            neighbors = list(self.graph.neighbors(stationId))
            if len(neighbors) != 2:
                raise Exception(f'Not exactly 2 neighbors connected to station "{ stationId }".')

            firstNeighbor = neighbors[0]
            edgeToFirstNeighbor = self.graph.get_edge_data(stationId, firstNeighbor)['edge']
            linesToFirstNeighbor = list(edgeToFirstNeighbor.getLines())

            secondNeighbor = neighbors[1]
            edgeToSecondNeighbor = self.graph.get_edge_data(stationId, secondNeighbor)['edge']
            linesToSecondNeighbor = list(edgeToSecondNeighbor.getLines())

            for lineToFirstNeighbor in linesToFirstNeighbor:
                if lineToFirstNeighbor not in linesToSecondNeighbor:
                    raise Exception(f'Lines along the two edges to station "{ stationId }" do not match.')

            for lineToSecondNeighbor in linesToSecondNeighbor:
                if lineToSecondNeighbor not in linesToFirstNeighbor:
                    raise Exception(f'Lines along the two edges to station "{ stationId }" do not match.')

            for line in linesToFirstNeighbor:
                self.removeConnection(stationId, firstNeighbor, line)
                self.removeConnection(stationId, secondNeighbor, line)

            stationData = self.graph.nodes[stationId]['station']
            self.removeStation(stationId)

            for line in linesToFirstNeighbor:
                self.addConnection(firstNeighbor, secondNeighbor, line)

            firstNeighborControlPoints = edgeToFirstNeighbor.getControlPoints(firstNeighbor)
            secondNeighborControlPoints = edgeToSecondNeighbor.getControlPoints(stationId)

            newEdge = self.graph.get_edge_data(firstNeighbor, secondNeighbor)['edge']
            newEdge.setControlPoints(
                firstNeighborControlPoints + [(stationData.x, stationData.y)] + secondNeighborControlPoints,
                firstNeighbor
            )

        return

    def addConnection(self, s1, s2, lineId, creator=None, information=None, controlPoints=[]):
        '''
        Adds a connection between to stations to the network. Adds an edge to the graph or adds the line to the existing edge of the graph.
        After adding an edge the iterator over the stations of a line is recalculated.

        :s1: First station
        :s2: Second station
        :lineId: The identifier of the line. (can also be list of IDs)
        :creator: Name of creator
        :information: Information added to the new edge, in case one is created
        :return: None
        '''

        if not isinstance(lineId, list):
            lineId = [lineId]

        for line in lineId:
            if line in self.lineGraphs:
                G = self.lineGraphs[line]
                G.add_edge(s1, s2)
            else:
                G = nx.Graph()
                G.add_edge(s1, s2)
                self.lineGraphs[line] = G
            self.recalculateLinearLine(line)

        if self.graph.has_edge(s1, s2):
            edge = self.graph.get_edge_data(s1, s2)['edge']
            for line in lineId:
                edge.addLine(line)
        else:
            edge = Edge(s1, s2, lineId, creator, information, controlPoints=controlPoints)
            self.graph.add_edge(s1, s2, edge=edge)

    def removeConnection(self, s1, s2, lineId):
        '''
        Removes a connection between to stations from the network. Removes the line from an existing edge of the graph and (if it is the last) also removes an edge from the graph.
        After removing an edge the iterator over the stations of a line is recalculated.

        :s1: First station
        :s2: Second station
        :lineId: The identifier of the line.

        :return: None
        '''
        if lineId in self.lineGraphs:
            G = self.lineGraphs[lineId]
            G.remove_edge(s1, s2)

        self.recalculateLinearLine(lineId)

        if self.graph.has_edge(s1, s2):
            edge = self.graph.get_edge_data(s1, s2)['edge']
            edge.removeLine(lineId)
            if len(edge.getLines()) == 0:
                self.graph.remove_edge(s1, s2)

    def convertIntersectionToStation(self, s1a, s2a, s1b, s2b, id, creatorName=None, information=None):
        if self.graph.has_edge(s1a, s2a) and self.graph.has_edge(s1b, s2b):
            station1aNode = self.graph.nodes[s1a]
            station2aNode = self.graph.nodes[s2a]
            edgeA = LineString([
                (station1aNode['station'].x, station1aNode['station'].y),
                (station2aNode['station'].x, station2aNode['station'].y)
            ])
            if list(edgeA.coords)[0] == list(edgeA.coords)[1]:
                return

            station1bNode = self.graph.nodes[s1b]
            station2bNode = self.graph.nodes[s2b]
            edgeB = LineString([
                (station1bNode['station'].x, station1bNode['station'].y),
                (station2bNode['station'].x, station2bNode['station'].y)
            ])
            if list(edgeB.coords)[0] == list(edgeB.coords)[1]:
                return

            if edgeA.intersects(edgeB):
                linesOfEdgeA = []
                linesOfEdgeB = []
                intersection = edgeA.intersection(edgeB)
                print(f"edgeA: {edgeA}, edgeB: {edgeB} intersect: {edgeA.intersects(edgeB)}")
                print(f"intersection: {intersection}")

                edgeAData = self.graph.get_edge_data(s1a, s2a)['edge']
                edgeBData = self.graph.get_edge_data(s1b, s2b)['edge']

                if len(edgeAData.controlPoints) > 0:
                    raise Exception(f'Edge "{ s1a } - { s2a }" must not contain control points. These would be removed when inserting stations.')

                if len(edgeBData.controlPoints) > 0:
                    raise Exception(f'Edge "{ s1b } - { s2b }" must not contain control points. These would be removed when inserting stations.')

                for line in edgeAData.getLines():
                    linesOfEdgeA.append(line)
                    self.removeConnection(s1a, s2a, line)

                for line in edgeBData.getLines():
                    linesOfEdgeB.append(line)
                    self.removeConnection(s1b, s2b, line)

                station = Station(id, id, intersection.x, intersection.y, creatorName, information)
                self.addStation(station)

                for line in linesOfEdgeA:
                    self.addConnection(s1a, station.id, line)
                    self.addConnection(s2a, station.id, line)

                self.graph.get_edge_data(s1a, station.id)['edge'].copyContextualInformationFrom(edgeAData)
                self.graph.get_edge_data(s2a, station.id)['edge'].copyContextualInformationFrom(edgeAData)

                for line in linesOfEdgeB:
                    self.addConnection(s1b, station.id, line)
                    self.addConnection(s2b, station.id, line)

                self.graph.get_edge_data(s1b, station.id)['edge'].copyContextualInformationFrom(edgeBData)
                self.graph.get_edge_data(s2b, station.id)['edge'].copyContextualInformationFrom(edgeBData)

    def getStation(self, u):
        return self.graph.nodes[u]['station']

    def getStationPosition(self, u):
        return self.graph.nodes[u]['station'].x, self.graph.nodes[u]['station'].y

    def getStations(self):
        return [data["station"] for _, data in self.graph.nodes(data=True)]

    def getEdge(self, s1, s2):
        '''
        Convienient getter method for edge objects.

        :s1: source
        :s2: target

        :return: Edge-object of (s1, s2)
        '''
        return self.graph[s1][s2]["edge"]

    def getEdges(self, s=None, data=False):
        '''
        Convenient getter method for edges adjacent to one station in the network.

        :s: station

        :return: (key, value) iterator over all edges
        '''
        return self.graph.edges(s, data)

    def recalculateLinearLine(self, line):
        '''
        Helper function to recalculate the linear order of stations along a line. Only works if all stations of a line are connected by an edge and the degree in regards to the line is <= 2 for all stations.

        :return: None
        '''
        G = self.lineGraphs[line]
        l = sorted(G.degree, key=lambda x: x[1], reverse=False)
        last = None
        for entry in l:
            if entry[1] > 0:
                last = entry[0]
                break
        if last is None:
            last = l[0][0]
        linearStations = [last]

        noneFound = False
        while not noneFound:
            noneFound = True

            for next in G.neighbors(last):
                if next not in linearStations:
                    linearStations.append(next)

                    last = next
                    noneFound = False

                    break

        self.lines[line] = linearStations

    def getMapDimensions(self, all_snapshots=False):
        maxx, maxy, minx, miny = -math.inf, -math.inf, math.inf, math.inf
        for station in self.getStations():
            minx = min(station.x, minx)
            miny = min(station.y, miny)
            maxx = min(station.x, maxx)
            maxy = min(station.y, maxy)
        return minx, miny, maxx, maxy


    def snapshot(self, name=''):
        '''
        Create a snapshot of the current map and store it in the map itself. 

        :param

        name: Identifier name of the snapshot
        '''
        existing_snapshots = self.snapshots
        self.snapshots = []
        copied_map = MapSnapshot(copy.deepcopy(self), name)
        self.snapshots = existing_snapshots
        self.snapshots.append(copied_map)

    def getInducedSubmap(self, stationIds):
        '''
        Does not use networkx' method in order to take advantage of our custom line handling
        :param stationIds: stations of the induced submap
        :return: the induced submap
        '''
        inducedSubmap = Map()
        for stationId in stationIds:
            if stationId in self.graph.nodes:
                inducedSubmap.addStation(copy.deepcopy(self.graph.nodes[stationId]['station']))

        for i in range(0, len(stationIds)):
            for j in range(i + 1, len(stationIds)):
                stationAId = stationIds[i]
                stationBId = stationIds[j]
                if stationAId != stationBId:
                    if self.graph.has_edge(stationAId, stationBId):
                        edge = copy.deepcopy(self.graph.get_edge_data(stationAId, stationBId)['edge'])
                        inducedSubmap.addConnection(stationAId, stationBId, edge.getLines())
                        copiedEdge = inducedSubmap.graph.get_edge_data(stationAId, stationBId)['edge']
                        copiedEdge.setControlPoints(edge.getControlPoints(stationAId), stationAId)
                        copiedEdge.copyContextualInformationFrom(edge)

        inducedSubmap.lineColors = copy.deepcopy(self.lineColors)

        return inducedSubmap

    def paste(self, other):
        '''
        Pastes a copy of another map into this map
        :param other: the other map

        In case a station id already exists, the station will not be added again, however new edges might be attached to the existing station
        '''
        for station_id in other.graph.nodes:
            if station_id not in self.graph.nodes:
                self.addStation(copy.deepcopy(other.graph.nodes[station_id]['station']))

        for u, v, data in other.graph.edges(data=True):
            edge_data = copy.deepcopy(data['edge'])
            self.addConnection(u, v, edge_data.getLines())
            copied_edge = self.graph.get_edge_data(u, v)['edge']
            copied_edge.setControlPoints(edge_data.getControlPoints(u), u)
            copied_edge.copyContextualInformationFrom(edge_data)

        for line in other.lineColors:
            if line not in self.lineColors:
                self.lineColors[line] = copy.deepcopy(other.lineColors[line])

    def reset_control_points(self):
        for id1, _, data in self.getEdges(data=True):
            data["edge"].setControlPoints([], id1)

    def set_station_positions(self, orig_positions):
        for s_id, data in self.graph.nodes(data=True):
            data["station"].x, data["station"].y = orig_positions[s_id]

    def get_station_positions(self):
        # new method 
        station_positions=[]
        for _ , data in self.graph.nodes(data=True):
            station_positions.append((data["station"].x,data["station"].y))
        return station_positions
            


class Station:
    '''
    Wrapper class for attributes regarding a station.
    '''


    def __init__(self, id, label=None, x=None, y=None, author=None, contextualInformation=None):
        '''
        Create a station. ID is required and all other attributes are optional. the ID is used as label if not specified.

        :id: ID of the station
        :label: Label of the station
        :x: x-coordinate of the station
        :y: y-coordinate of the station
        :author: Name of author --- shortcut to set contextual information on initialization
        :contextualInformation: Contextual information added by an author --- shortcut to set contextual information on initialization
        '''
        self.id = id
        self.x = x
        self.y = y

        self.contextualInformation = ContextualInformation()
        if author is not None:
            self.contextualInformation.setInformation(author, contextualInformation)

        if label:
            self.label = label
        else:
            self.label = self.id


    def hasContextualInformation(self, author):
        '''
        :return: bool if contextual information of author exists
        '''
        return self.contextualInformation.hasInformation(author)

    def getContextualInformation(self, author):
        '''
        :return: get contextual information of author
        '''
        return self.contextualInformation.getInformation(author)

    def setContextualInformation(self, author, information):
        '''
        set contextual information for author
        '''
        self.contextualInformation.setInformation(author, information)

    def removeContextualInformation(self, author):
        '''
        remove contextual information of author
        '''
        self.contextualInformation.removeInformation(author)

    def copyContextualInformationFrom(self, station, override=True):
        '''
        copies the contextual information of other station
        :override: if true, existing information of the same author will be overridden
        '''
        authorsOfOtherStation = station.contextualInformation.getAuthors()
        for author in authorsOfOtherStation:
            if self.hasContextualInformation(author):
                if override:
                    self.setContextualInformation(author, copy.deepcopy(station.getContextualInformation(author)))
            else:
                self.setContextualInformation(author, copy.deepcopy(station.getContextualInformation(author)))


class Edge:
    '''
    Wrapper class for an edge in the map.
    '''
    source = None
    target = None
    sourcePort = []
    targetPort = []
    controlPoints = []

    def __init__(self, source, target, line, author=None, contextualInformation=None, controlPoints=[]):
        '''
        Create an edge that wraps edge attributes. Source and target order is not important.

        :source: The ID of the source station
        :target: The ID of the target station
        :line: The ID of the line that uses the respective edge (can also be list of IDs)
        :author: Name of author --- shortcut to set contextual information on initialization
        :contextualInformation: Contextual information added by an author --- shortcut to set contextual information on initialization
        '''
        self.source = source
        self.target = target
        if isinstance(line, list):
            self.sourcePort = line.copy()
            self.targetPort = line.copy()
        else:
            self.sourcePort = [line]
            self.targetPort = [line]

        self.contextualInformation = ContextualInformation()
        if author is not None:
            self.contextualInformation.setInformation(author, contextualInformation)
        self.controlPoints = controlPoints
        return


    # def getEdgeGeometry(self, start=None):
    #     '''
    #     Returns the actual geometry of the routed edge at the moment (source -> control points -> target)
    #
    #     :return: [src, (ctrl_points)*, tgt] list of geometry points of the routed edge
    #     '''
    #     if start == self.source or start == None:
    #         return [self.source] + self.controlPoints + [self.target]
    #     elif start == self.target:
    #         return [self.target] + self.controlPoints[::-1] + [self.source]
    #     else:
    #         print(f"provided start {start} is neither source nor target of edge {self.source}-{self.target}")

    def getControlPoints(self, start=None):
        '''
        Returns the control points of the routed edge at the moment (first control point is nearest to end specified in 'start')

        :return: controlPoints list of geometry points of the routed edge
        '''
        if start == self.source or start == None:
            return self.controlPoints
        elif start == self.target:
            return self.controlPoints[::-1]
        else:
            print(f"provided start {start} is neither source nor target of edge {self.source}-{self.target}")

    def setControlPoints(self, controlPoints, start=None):
        '''
        Sets the control points of the edge (first control point is nearest to end specified in 'start')

        :return: None
        '''
        if start == self.source or start == None:
            self.controlPoints = controlPoints
        elif start == self.target:
            self.controlPoints = controlPoints[::-1]
        else:
            print(f"provided start {start} is neither source nor target of edge {self.source}-{self.target}")

    def addLine(self, line):
        '''
        Adds a line that uses an edge in the network. 

        :line: the identifier of the line that is added to the edge.
        '''
        self.sourcePort.append(line)
        self.targetPort.append(line)

    def removeLine(self, line):
        '''
        Removes a line that uses an edge in the network.

        :line: the identifier of the line that is removed from the edge.
        '''
        self.sourcePort.remove(line)
        self.targetPort.remove(line)

    def getPort(self, station):
        '''
        Get the ports an edge
        '''
        if station == self.source:
            return self.sourcePort
        if station == self.target:
            return self.targetPort

    def getLines(self):
        '''
        Get the list of lines on this edge (uses source port, order is therefore not necessarily useful)
        '''
        return self.sourcePort

    def hasContextualInformation(self, author):
        '''
        :return: bool if contextual information of author exists
        '''
        return self.contextualInformation.hasInformation(author)

    def getContextualInformation(self, author):
        '''
        :return: get contextual information of author
        '''
        return self.contextualInformation.getInformation(author)

    def setContextualInformation(self, author, information):
        '''
        set contextual information for author
        '''
        self.contextualInformation.setInformation(author, information)

    def removeContextualInformation(self, author):
        '''
        remove contextual information of author
        '''
        self.contextualInformation.removeInformation(author)

    def copyContextualInformationFrom(self, edge, override=True):
        '''
        copies the contextual information of other edge
        :override: if true, existing information of the same author will be overridden
        '''
        authorsOfOtherEdge = edge.contextualInformation.getAuthors()
        for author in authorsOfOtherEdge:
            if self.hasContextualInformation(author):
                if override:
                    self.setContextualInformation(author, copy.deepcopy(edge.getContextualInformation(author)))
            else:
                self.setContextualInformation(author, copy.deepcopy(edge.getContextualInformation(author)))

