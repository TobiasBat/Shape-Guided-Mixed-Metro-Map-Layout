import math

from ac_metro.utility.abstractUtility import AbstractUtility


class GridPlacementUtility(AbstractUtility):
    '''
    Compute the best candidate set (might be singular vertex) for a station in an embedded graph
    '''

    @staticmethod
    def execute(map, options=None):

        # the grid in which we find the best candidates for each station of the map
        #   - can be irregular
        #   - must have x and y coordinates as properties of the nodes
        grid = options["grid"]

        # prefer grid nodes with this attribute set to true
        if "prefer_attribute" in options:
            prefer_att = options["prefer_attribute"]
            prefer_factor = options["prefer_factor"]
            print(f"preferring {prefer_att}")

        # argument with varying meaning depending on the method
        k = options["k"]

        # size of a grid cell
        grid_cell_size = options["cell_size"]

        # possibilities are:
        #   knearest:   set of k closest neighbors (set k to 1 for exactest) [naive implementation]
        #   circle:     set of all grid positions inside a circle of radius k
        #   diamond:    set of all grid positions with L_1 distance <= k
        #   rectangle:  set of all grid positions with L_{inf} distance <= k
        method = options["method"]

        # Ignore all vertices, with a degree smaller than this
        min_deg = 0
        if "min_degree" in options:
            min_deg = options["min_degree"]

        # Ignore all vertices, with this property set to true
        ignore = ""
        if "ignore" in options:
            ignore = options["ignore"]

        # if the "custom" method is chosen, the passed distance function is used to determine the distance
        # function must have the signature dist(x1, y1, x2, y2) and return a number
        if "dist_func" in options:
            distance_function = options["dist_func"]

        candidates = {}
        distances = {}

        if method == "knearest":
            for i, (node, data) in enumerate(map.graph.nodes(data=True)):
                station = data["station"]
                nodes_distances = []
                for j, (g_node, g_data) in enumerate(grid.nodes(data=True)):

                    #Skipping vertices with too few adjacencies
                    if grid.degree[g_node]<min_deg:
                        continue
                    if ignore != "" and ignore in g_data and g_data[ignore]:
                        # print(f"skipping {g_node}")
                        continue

                    nodes_distances.append(
                        (g_node,
                         math.sqrt(
                             math.pow(station.x - g_data['x'], 2) +
                             math.pow(station.y - g_data['y'], 2)
                         ))
                    )
                    distances[node] = sorted(nodes_distances, key=lambda x: x[1])
            for node in distances:
                candidates[node] = [i[0] for i in distances[node][:k]]
        else:
            for i, (node, data) in enumerate(map.graph.nodes(data=True)):
                station = data["station"]
                nodes_candidates = []
                for j, (g_node, g_data) in enumerate(grid.nodes(data=True)):
                    
                    #Skipping vertices with too few adjacencies
                    if grid.degree[g_node]<min_deg:
                        continue
                    if ignore != "" and ignore in g_data and g_data[ignore]:
                        # print(f"skipping {g_node}")
                        continue

                    if method == "diamond":
                        distance = (station.x - g_data['x']) + (station.y - g_data['y'])
                    elif method == "circle":
                        distance = math.sqrt(
                            math.pow(station.x - g_data['x'], 2) +
                            math.pow(station.y - g_data['y'], 2)
                        )
                    elif method == "square":
                        distance = max(station.x - g_data['x']), (station.y - g_data['y'])
                    elif method == "custom":
                        distance = distance_function(station.x, station.y, g_data['x'], g_data['y'])
                    else:
                        print("No valid method has been provided to GridPlacementUtility")
                    if distance <= k * grid_cell_size:
                        # print(f"info {grid.nodes()[g_node]}")
                        # input()
                        if "prefer_attribute" in options and grid.nodes()[g_node][prefer_att]:
                            distance = distance/prefer_factor # Here we pretende like this node is closer (after checking the actual distance)
                            # print(f"\tshortened distance for {g_node}")
                        # print(f"\tappending {g_node}")
                        nodes_candidates.append((g_node, distance))

                # we return either just an ordered list of candidates as value for the key, which is an ID of a station
                # or we return an ordered list of tuples, the first being again the position and the second being the value they are sorted by (e.g. distance to optimum)
                if "distances" in options and options["distances"]:
                    candidates[node] = nodes_candidates
                else:
                    candidates[node] = [i[0] for i in nodes_candidates]
        return candidates
