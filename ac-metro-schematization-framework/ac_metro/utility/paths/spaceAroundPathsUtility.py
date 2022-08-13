from shapely.geometry import LineString, Point

from ac_metro.utility.abstractUtility import AbstractUtility
from ac_metro.utility.paths.getPathsUtility import GetPathsUtility
from wktplot import WKTPlot


def get_path_geometry(Map, path):
    stations = [path.getStartId()] + [station.id for station in path.getInnerStations()] + [path.getEndId()]
    geometry = [Map.getStationPosition(path.getStartId())]

    for i in range(1, len(stations)):
        edge_object = Map.getEdge(stations[i-1], stations[i])
        # print(f"edge {(stations[i-1], stations[i])} has cpoints: {edge_object.getControlPoints(start=stations[i - 1])}")
        geometry += edge_object.getControlPoints(start=stations[i-1])
        geometry.append(Map.getStationPosition(stations[i]))
        edge_object.setControlPoints([])

    # edges = path.getInnerEdges()
    # print(edges)
    # IDs = path.getInnerIDs()
    # last_station = path.getStartId()
    # for i, edge in enumerate(edges):
    #     edge_object = edge["edge"]
    #     print(f"edge {(edge_object.source, edge_object.target)} with id {IDs[i]}")
    #     geometry += edge_object.getControlPoints(start=last_station)
    #     last_station = IDs[i]
    #     edge_object.setControlPoints([])
    # geometry += [Map.getStationPosition(path.getEndId())]
    return geometry


class SpaceAroundPathsUtility(AbstractUtility):
    '''
    Repositions inner vertices of paths s.t. there is equal space around them
    '''

    @staticmethod
    def UTIL_NAME():
        return 'SpaceAroundPathsUtility'

    @staticmethod
    def execute(Map, options={}):
        path_information = GetPathsUtility.execute(Map)
        for path in path_information.getPaths():
            path_geom = get_path_geometry(Map, path)
            # print(f"And its geometry is {path_geom}")
            line_string = LineString(path_geom)
            # Add shapes to the plot
            control_point_positions = list(map(lambda controlpoint: (line_string.project(Point(controlpoint), True), controlpoint), path_geom))
            control_point_positions.sort(key=(lambda x: x[0]))
            print(control_point_positions)
            # path_length = line_string.length
            steps = float(len(path.getInnerStations())+1)

            stations = [path.getStartId()] + [station.id for station in path.getInnerStations()] + [path.getEndId()]
            control_point_pointer = 0
            last_id = path.getStartId()
            last_fraction = 0

            for i in range(1, len(stations)):
            # for i, station in enumerate(path.getInnerStations()):
                station = Map.getStation(stations[i])
                new_id = station.id
                fraction = float(i) / steps
                new_pos = line_string.interpolate(fraction, True)
                station.x = new_pos.x
                station.y = new_pos.y
                new_control_points = []
                for cpoint in control_point_positions:
                    if last_fraction < cpoint[0] < fraction:
                        new_control_points.append(cpoint[1])
                # while control_point_positions[control_point_pointer][0] < fraction:
                #     new_control_points.append(control_point_positions[control_point_pointer][1])
                #     control_point_pointer += 1
                Map.getEdge(last_id, new_id).setControlPoints(new_control_points, start=last_id)
                last_id = new_id
                last_fraction = fraction




    @staticmethod
    def execute_old(Map, options={}):

        path_information = GetPathsUtility.execute(Map)
        for path in path_information.getPaths():
            print(f"path: {path}")
            inner_stations = path.getInnerStations()
            first_station = Map.graph.nodes[path.getStartId()]['station']
            last_station = Map.graph.nodes[path.getEndId()]['station']
            print(f"first: {first_station.id}")
            print(f"crossing: {[st.id for st in inner_stations]}")
            print(f"last: {last_station.id}")
            stations = [first_station] + inner_stations + [last_station]
            stations_to_be_replaced_with_controlpoints = []
            for i in range(1, len(stations) - 1):

                station = stations[i]
                previous_station = stations[i - 1]
                next_station = stations[i + 1]
                line = LineString([
                    (previous_station.x, previous_station.y),
                    (next_station.x, next_station.y)
                ])
                if line.intersects(Point(station.x, station.y)):
                    Map.removeStationKeepLines(station.id)
                else:
                    stations_to_be_replaced_with_controlpoints.append(station)

            for station in stations_to_be_replaced_with_controlpoints:
                Map.convertStationToControlPoint(station.id)

            Map.addStationsBetween(
                path.getStartId(),
                path.getEndId(),
                len(inner_stations),
                list(map(lambda station: station.id, inner_stations)),
                None,
                None,
                None)

