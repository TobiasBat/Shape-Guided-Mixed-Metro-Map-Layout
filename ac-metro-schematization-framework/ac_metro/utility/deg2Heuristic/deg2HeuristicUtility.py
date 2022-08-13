import math

from ac_metro.map.map import Station
from ac_metro.utility.abstractUtility import AbstractUtility
from ac_metro.utility.paths.getPathsUtility import GetPathsUtility
from copy import deepcopy


class Deg2HeuristicInformation:
    def __init__(self, numberOfReplacedEdges, removedEdges=None):
        self.numberOfReplacedEdges = numberOfReplacedEdges
        self.removedEdges = removedEdges


class Deg2HeuristicUtility(AbstractUtility):
    '''
    Adjust the amount of nodes in all paths
    '''

    @staticmethod
    def UTIL_NAME():
        return 'Deg2HeuristicUtility'

    @staticmethod
    def execute(Map, options={}):
        count = 0

        number_of_nodes_per_path = 2
        if 'number_of_nodes_per_path' in options and options['number_of_nodes_per_path'] >= 0:
            if options['number_of_nodes_per_path'] == 0:
                print(f'Setting "number_of_nodes_per_path" of {Deg2HeuristicUtility.UTIL_NAME()} could result in loss of edges.')
            number_of_nodes_per_path = options['number_of_nodes_per_path']

        number_of_nodes_per_path_to_terminal_stations = 1
        if 'number_of_nodes_per_path_to_terminal_stations' in options and options['number_of_nodes_per_path_to_terminal_stations'] >= 0:
            if options['number_of_nodes_per_path_to_terminal_stations'] == 0:
                print(f'Setting "number_of_nodes_per_path_to_terminal_stations" of {Deg2HeuristicUtility.UTIL_NAME()} could result in loss of edges.')
            number_of_nodes_per_path_to_terminal_stations = options['number_of_nodes_per_path_to_terminal_stations']

        do_not_add_new_stations = False
        if 'do_not_add_new_stations' in options:
            do_not_add_new_stations = options['do_not_add_new_stations']

        stations_that_have_to_remain = []
        if 'stations_that_have_to_remain' in options:
            stations_that_have_to_remain = options['stations_that_have_to_remain']

        only_process_paths = None
        if 'only_process_paths' in options:
            only_process_paths = options['only_process_paths']

        added_helper_single_degree_stations = []
        count_added_helper_single_degree_stations = 0
        for station_that_has_to_remain in stations_that_have_to_remain:
            id_added_helper_single_degree_station_1 = f'deg2heu_internal_single_degree_helper_{count_added_helper_single_degree_stations}'
            id_added_helper_single_degree_station_2 = f'deg2heu_internal_single_degree_helper_{count_added_helper_single_degree_stations + 1}'
            Map.addStation(Station(id_added_helper_single_degree_station_1, None, 0, 0, Deg2HeuristicUtility.UTIL_NAME()))
            Map.addStation(Station(id_added_helper_single_degree_station_2, None, 0, 0, Deg2HeuristicUtility.UTIL_NAME()))
            Map.addConnection(station_that_has_to_remain, id_added_helper_single_degree_station_1, f'deg2heu_internal_single_stop_helper_line_{count_added_helper_single_degree_stations}', Deg2HeuristicUtility.UTIL_NAME())
            Map.addConnection(station_that_has_to_remain, id_added_helper_single_degree_station_2, f'deg2heu_internal_single_stop_helper_line_{count_added_helper_single_degree_stations + 1}', Deg2HeuristicUtility.UTIL_NAME())
            count_added_helper_single_degree_stations += 2
            added_helper_single_degree_stations.append(id_added_helper_single_degree_station_1)
            added_helper_single_degree_stations.append(id_added_helper_single_degree_station_2)

        path_information = GetPathsUtility.execute(Map)
        for path in path_information.getPaths():

            if only_process_paths is not None:
                included = False
                for only_process_path in only_process_paths:
                    if (path.getStartId() == only_process_path[0] and path.getEndId() == only_process_path[1]) \
                            or (path.getStartId() == only_process_path[1] and path.getEndId() == only_process_path[0]):
                        included = True
                        break
                if not included:
                    continue

            if path.getStartId() in added_helper_single_degree_stations or path.getEndId() in added_helper_single_degree_stations:
                continue

            if path.getStartId() in stations_that_have_to_remain \
                    and path.getEndId() in stations_that_have_to_remain \
                    and len(path.getInnerStations()) == 0:
                continue

            number_of_nodes = number_of_nodes_per_path

            numberOfNeighborsOfPathStartNode = len(list(Map.graph.neighbors(path.getStartId())))
            numberOfNeighborsOfPathEndNode = len(list(Map.graph.neighbors(path.getEndId())))
            if numberOfNeighborsOfPathStartNode <= 1 or numberOfNeighborsOfPathEndNode <= 1:
                number_of_nodes = number_of_nodes_per_path_to_terminal_stations

            if len(path.getInnerStations()) < number_of_nodes and not do_not_add_new_stations:
                if len(path.getInnerStations()) == 0:
                    ids_of_inserted_stations = list(map(lambda x: f"deg2heu_{x}", list(range(count, count + number_of_nodes))))
                    Map.addStationsBetween(
                        path.getStartId(),
                        path.getEndId(),
                        number_of_nodes,
                        ids_of_inserted_stations,
                        Deg2HeuristicUtility.UTIL_NAME(),
                        None,
                        ([Deg2HeuristicInformation(0, None)] * number_of_nodes) + [Deg2HeuristicInformation(1, None)])
                    count += number_of_nodes

                elif len(path.getInnerStations()) > 0:
                    difference = number_of_nodes - len(path.getInnerStations())

                    firstHalf = math.floor(difference / 2)
                    if numberOfNeighborsOfPathEndNode <= 1: # if terminal station is on "end" side, move more new stations toward "start" side
                        firstHalf = math.ceil(difference / 2)
                    ids_of_inserted_stations = list(map(lambda x: f"deg2heu_{x}", list(range(count, count + firstHalf))))
                    Map.addStationsBetween(
                        path.getStartId(),
                        path.getInnerStations()[0].id,
                        firstHalf,
                        ids_of_inserted_stations,
                        Deg2HeuristicUtility.UTIL_NAME(),
                        None,
                        ([Deg2HeuristicInformation(0, None)] * firstHalf) + [Deg2HeuristicInformation(1, None)])
                    count += firstHalf

                    secondHalf = difference - firstHalf
                    ids_of_inserted_stations = list(map(lambda x: f"deg2heu_{x}", list(range(count, count + secondHalf))))
                    Map.addStationsBetween(
                        path.getInnerStations()[-1].id,
                        path.getEndId(),
                        secondHalf,
                        ids_of_inserted_stations,
                        Deg2HeuristicUtility.UTIL_NAME(),
                        None,
                        ([Deg2HeuristicInformation(0, None)] * secondHalf) + [Deg2HeuristicInformation(1, None)])
                    count += secondHalf

            elif len(path.getInnerStations()) > number_of_nodes:
                firstHalf = math.floor(number_of_nodes / 2)
                secondHalf = math.ceil(number_of_nodes / 2)
                if numberOfNeighborsOfPathEndNode <= 1: # if terminal station is on "end" side, keep more stations at "start" side
                    firstHalf = math.ceil(number_of_nodes / 2)
                    secondHalf = math.floor(number_of_nodes / 2)
                removedLines = []
                removedEdges = deepcopy(path.getInnerEdges()[firstHalf:(len(path.getInnerEdges()) - secondHalf)])
                for i in range(firstHalf, len(path.getInnerStations()) - secondHalf):
                    stationId = path.getInnerStations()[i].id
                    neighbors = list(Map.graph.neighbors(stationId))

                    for neighbor in neighbors:
                        edge = Map.graph.get_edge_data(stationId, neighbor)['edge']
                        lines = edge.getLines()
                        for line in list(lines):
                            if line not in removedLines:
                                removedLines.append(line)
                            Map.removeConnection(stationId, neighbor, line)

                    Map.removeStation(stationId)

                previousStationRemaining = None
                if firstHalf - 1 >= 0:
                    previousStationRemaining = path.getInnerStations()[firstHalf - 1].id
                else:
                    previousStationRemaining = path.getStartId()

                nextStationRemaining = None
                if len(path.getInnerStations()) - secondHalf < len(path.getInnerStations()):
                    nextStationRemaining = path.getInnerStations()[len(path.getInnerStations()) - secondHalf].id
                else:
                    nextStationRemaining = path.getEndId()

                Map.addConnection(previousStationRemaining, nextStationRemaining, removedLines, Deg2HeuristicUtility.UTIL_NAME(), Deg2HeuristicInformation(len(removedEdges), removedEdges))

        for added_helper_single_degree_station in added_helper_single_degree_stations:
            neighbors = list(Map.graph.neighbors(added_helper_single_degree_station))
            if len(neighbors) != 1:
                raise Exception(f'Error with helper single degree station {added_helper_single_degree_station}.')
            edge = Map.graph.get_edge_data(added_helper_single_degree_station, neighbors[0])['edge']
            for line in list(edge.getLines()):
                Map.removeConnection(added_helper_single_degree_station, neighbors[0], line)
            Map.removeStation(added_helper_single_degree_station)

        return

