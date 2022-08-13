# A helper class keeping an extended version of the grid including metanodes for modeling turn costs
import itertools
import math
import networkx as nx
import numpy as np

from ac_metro.utility.custom.octi.helper import dist, compute_x_axis_angle, compute_bend_angle, angle_cost, angle_cost_2, S_HOOK

DEBUG = False

class ExtendedGrid:

    def __init__(self, simpleGrid, bend_penalty) -> None:
        '''
        Given any graph, every node is replaced with meta nodes (cliques). Coordinates of port nodes are only for drawing the grid.
        Sets relevant corssings and port <-> parent mappings as well as next and previous port maps.
        '''
        self.sinkCost = 1000  # This should be a large constant, and technically it should be based on the input weights TODO: Change
        self.simpleGrid = simpleGrid
        self.graph = nx.Graph()
        self.starting_weights = {}  # A dict keeping track of initial values of edges (for resetting)
        self.port_mapping = {}  # This tells us by what port connection an edge will be replaced
        self.bend_penalty = bend_penalty

        # Sanity Check
        for v, data in self.simpleGrid.nodes(data=True):
            if "shape_part" not in data:
                print("Extended grid is created from a graph, in which not every node has  shape_part attribute.")

        # Adding original nodes of simple grid TODO: this is kinda ugly
        # D: could be max(list(self.graph.nodes))
        firs_usable_id = self.add_original_nodes()

        print(f"Extended grid created with {simpleGrid.number_of_nodes()} nodes in simple grid and first usable ID is {firs_usable_id}")

        # Calculating the smallest distance between two points in the simple grid
        smallest = self.find_shortest_distance(simpleGrid)

        # We use this smallest distance to draw port nodes legibly.
        # This is also functionally important for the port order, which is computed geometrically right now
        # TODO: This shouldn't be done geometrically
        movement = smallest

        # Adding ports
        self.add_port_nodes_and_connections(firs_usable_id, movement)

        # Adding crossing (based on crossing in simple grid)
        self.add_crossings()

        # saving next and previous ports (radial order) and parents
        self.port_next, self.port_prev, self.port_prnt, self.ports = {}, {}, {}, {}

        # Saving the previous next and parent/ports for all nodes
        self.mark_ports_and_parents()

        # Create connections between ports of the same node. This is where bend penalties are applied
        self.create_inter_port_connections()

    # Construction methods
    # HEREEEE START !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    def create_inter_port_connections(self):
        for u, data in self.simpleGrid.nodes(data=True):
            for port1, port2 in itertools.combinations(self.ports[u], 2):
                if port1 == port2:
                    print("CHANGE SOMETHING")
                angle = compute_bend_angle(
                    self.graph.nodes[port1]['x'], self.graph.nodes[port1]['y'],
                    self.graph.nodes[u]['x'], self.graph.nodes[u]['y'],
                    self.graph.nodes[port2]['x'], self.graph.nodes[port2]['y'])
                if data['shape_part']:
                    cost = 0
                else:
                    cost = angle_cost_2(angle) * self.bend_penalty
                # cost = angle_cost_2(angle) * self.bend_penalty
                # print(f"angle between {port1} and {port2} is {math.degrees(angle)} and costs {cost}")
                # input()

                self.graph.add_edge(port1, port2, weight=cost, crossing=[])

                self.starting_weights[frozenset((port1, port2))] = cost
                self.starting_weights[frozenset((port2, port1))] = cost

    def mark_ports_and_parents(self):
        for u in self.simpleGrid.nodes():
            nghbs = [n for n in self.graph.neighbors(u)]
            nghbs.sort(key=lambda x: compute_x_axis_angle(self.graph.nodes[u]['x'], self.graph.nodes[u]['y'],
                                                          self.graph.nodes[x]['x'], self.graph.nodes[x]['y']))

            # # Sanity Check
            # print(f"{u} has {len(nghbs)} neighbors:\n{nghbs}\nwhich are set as its ports")

            num_nghbs = len(nghbs)
            # print(f"when setting the parent {u} has {num_nghbs} ports, i.e., {nghbs}")
            self.ports[u] = nghbs
            for i in range(num_nghbs):
                self.port_next[nghbs[i]] = nghbs[(i + 1) % num_nghbs]
                self.port_prev[nghbs[i]] = nghbs[(i - 1) % num_nghbs]
                self.port_prnt[nghbs[i]] = u

            # Sanity Check
            for port in nghbs:
                if self.get_parent(port) != self.get_parent(self.port_next[port]) or self.get_parent(port) != self.get_parent(self.port_prev[port]):
                    print("Somethings weird")
                    input()
                # print(f"the port {port} has\nprevious: {self.port_prev[port]}\nnext: {self.port_next[port]}\nparent: {self.port_prnt[port]}")

    def add_crossings(self):
        for us, vs, data in self.simpleGrid.edges(data=True):
            u, v = self.port_mapping[frozenset((us, vs))]
            crossing_edges = data["crossing"]
            if len(crossing_edges) == 0:
                continue
            for edge in crossing_edges:
                u1, v1 = self.port_mapping[frozenset((edge[0], edge[1]))]
                self.graph[u][v]["crossing"].append((u1, v1))

    def find_shortest_distance(self, simpleGrid):
        smallest = math.inf
        for u, v, data in simpleGrid.edges(data=True):
            distance = dist(self.graph.nodes[u]['x'], self.graph.nodes[u]['y'], self.graph.nodes[v]['x'],
                            self.graph.nodes[v]['y'])
            # distance = math.sqrt( math.pow(self.graph[u]['x']-self.graph[v]['x'], 2) + math.pow(self.graph[u]['y']-self.graph[v]['y'], 2) )
            if distance < smallest:
                smallest = distance
        return smallest

    def add_original_nodes(self):
        # D: could use the function self.graph.add_nodes_from(self.simpleGrid.nodes(data=True))
        max_id = 0
        for u, data in self.simpleGrid.nodes(data=True):
            self.graph.add_node(u, x=data["x"], y=data["y"],
                                shape_part=data["shape_part"] if "shape_part" in data else False)
            if u > max_id:
                max_id = u
        return max_id + 1

    def add_port_nodes_and_connections(self, max_id, movement):
        for u, v, data in self.simpleGrid.edges(data=True):
            # weight of the edge in the simple grid
            starting_weight = data["weight"]
            # print(f"starting_weight of {(u, v)}: {starting_weight}")

            # Calculating the offset for a single port node
            direction = np.array([
                (self.graph.nodes[u]['x'] - self.graph.nodes[v]['x']),
                (self.graph.nodes[u]['y'] - self.graph.nodes[v]['y'])
            ])
            direction = direction / np.sqrt(np.sum(direction ** 2))
            # direction *= (1.0 * movement)
            direction *= (0.2 * movement)

            # Adding two ports on this edge (SHOULD SHAPE INFO BE ADDED?)
            self.graph.add_node(max_id, x=self.graph.nodes[u]['x'] - direction[0],
                                y=self.graph.nodes[u]['y'] - direction[1])
            self.graph.add_node(max_id + 1, x=self.graph.nodes[v]['x'] + direction[0],
                                y=self.graph.nodes[v]['y'] + direction[1])

            # Connecting node - port - port - node
            self.graph.add_edge(u, max_id, weight=self.sinkCost, crossing=[])
            self.graph.add_edge(max_id + 1, v, weight=self.sinkCost, crossing=[])
            self.graph.add_edge(max_id, max_id + 1, weight=starting_weight, crossing=[])

            if (u == 129 and max_id == 9897) or (v == 129 and max_id+1 == 9897):
                print(f"edge (129, 9897) exists")
                input()

            self.starting_weights[frozenset((u, max_id))] = self.sinkCost
            self.starting_weights[frozenset((max_id, u))] = self.sinkCost
            self.starting_weights[frozenset((max_id + 1, v))] = self.sinkCost
            self.starting_weights[frozenset((v, max_id + 1))] = self.sinkCost
            self.starting_weights[frozenset((max_id + 1, max_id))] = starting_weight
            self.starting_weights[frozenset((max_id, max_id + 1))] = starting_weight
            # print(f"{u}, {max_id}, {max_id+1}, {v}")

            self.port_mapping[frozenset((u, v))] = (max_id, max_id + 1)

            max_id += 2
        print(f"added nodes for extended grid up to {max_id-1}")
    # End construction methods


    # Accessible methods
    def coord(self, node):
        '''
        returns the coordinates of a grid node
        '''
        if self.is_port(node):
            return self.graph.nodes[self.get_parent(node)]["x"], self.graph.nodes[self.get_parent(node)]["y"]
        else:
            return self.graph.nodes[node]["x"], self.graph.nodes[node]["y"]

    def is_port(self, ID):
        '''
        Checks if node is a port
        '''
        return ID not in self.ports

    def remove_ports(self, path):
        '''
        Given a path as a list of nodes, removes all ports from the list.
        '''
        simple_path = []
        for node in path:
            if not self.is_port(node):
                simple_path.append(node)
        return simple_path

    def reset_edge_weights(self, edge=None):
        '''
        Resets the weight of every edge back to its inital cost
        '''
        print(f"Resetting edge weights")
        # print(self.starting_weights.keys())
        if edge == None:
            for u, v in self.graph.edges():
                try:
                    self.graph[u][v]['weight'] = self.starting_weights[frozenset((u, v))]
                except KeyError:
                    print(f"looking for {(u, v)}, which is {'' if self.graph.has_edge(u, v) else 'not'} in self.graph.")
                    input()
        else:
            self.graph[edge[0]][edge[1]]['weight'] = self.starting_weights[frozenset((edge[0], edge[1]))]

    def close_path(self, path):
        '''
        Closes node or parent (if node is a port) and all port to port edges, by setting them to math.inf. Argument is a list of nodes.
        '''
        # print(f"closing {path}")
        for i in range(len(path)):
            if not self.is_port(path[i]):
                self.close_node(path[i])
            else:
                try:
                    self.close_node(self.port_prnt[path[i]])
                except KeyError:
                    print(f"{path[i]} has no parent")
                    input()

            # for edge in self.graph.edges(path_node):
            #     self.close_edge(edge[0], edge[1])
            #     # self.graph[edge[0]][edge[1]]["weight"] = math.inf

            if i > 0:
                if self.is_port(path[i]) and self.is_port(path[i-1]):
                    self.close_edge(path[i-1], path[i])
                for cross_edge in self.graph[path[i - 1]][path[i]]["crossing"]:
                    self.close_edge(cross_edge[0], cross_edge[1])
                    # self.graph[cross_edge[0]][cross_edge[1]]["weight"] = math.inf

    def close_node(self, node):
        '''
        Closes all sink nodes and connections between ports of a node
        '''
        if self.is_port(node):
            print(f"{node} is a port node. Close the parent instead!")
        else:
            for port in self.ports[node]:
                self.close_edge(node, port)
                # self.graph[node][port]["weight"] = math.inf
            for port1, port2 in itertools.combinations(self.ports[node], 2):
                self.close_edge(port1, port2)
                # self.graph[port1][port2]["weight"] = math.inf

    def open_node(self, node):
        '''
        Opens all sink nodes and connections between ports of a node
        '''
        if not self.is_port(node):
            for port in self.ports[node]:
                self.open_edge(node, port)
            for port1, port2 in itertools.combinations(self.ports[node], 2):
                self.open_edge(port1, port2)

    def open_all_ports(self, pos):
        for port in self.ports[pos]:
            self.open_edge(pos, port)


    def open_ports(self, src_pos, next_edge_pos, prev_edge_pos, next_off, prev_off):
        '''
        Opens all sink edges from one edge to another minus offsets on both sides
        '''
        if DEBUG:
            print(f"open ports called with")
            print(f"srcpos: {src_pos} with ports {self.ports[src_pos]}")
            print(f"next_edge_pos: {next_edge_pos}")
            print(f"prev_edge_pos: {prev_edge_pos}")
            print(f"next_off: {next_off}")
            print(f"prev_off: {prev_off}")

        # print(f"one step from {prev_edge_pos} to {self.port_next[prev_edge_pos]} for start")
        start_port = self.port_prev[prev_edge_pos]
        for i in range(prev_off):
            start_port = self.port_prev[start_port]

        # print(f"one step from {prev_edge_pos} to {self.port_prev[prev_edge_pos]} for end")
        end_port = self.port_next[next_edge_pos]
        for i in range(next_off):
            end_port = self.port_next[end_port]

        if DEBUG:
            print(f"before opening ports start: {start_port} and end: {end_port}")
        while start_port != end_port:
            if DEBUG:
                print(f"start port: {start_port} and end port: {end_port}")
            self.open_edge(src_pos, start_port)
            start_port = self.port_prev[start_port]
        self.open_edge(src_pos, start_port)

    def open_path(self, path):
        for i in range(0, len(path)-1):
            if not self.is_port(path[i]):
                self.open_node(path[i])
            else:
                self.open_node(self.port_prnt[path[i]])

            if i > 0:
                if self.is_port(path[i]) and self.is_port(path[i-1]):
                    self.close_edge(path[i-1], path[i])
                for cross_edge in self.graph[path[i - 1]][path[i]]["crossing"]:
                    self.open_edge(cross_edge[0], cross_edge[1])

    def open_edge(self, u, v, offset=0):
        '''
        Opens an edge, by resetting the weight to its original
        '''
        # print(f"opening edge ({u}-{v})")
        self.graph[u][v]['weight'] = self.starting_weights[frozenset((u, v))] + offset
        for cross_edge in self.graph[u][v]["crossing"]:
            a, b = cross_edge[0], cross_edge[1]
            self.graph[a][b]['weight'] = self.starting_weights[frozenset((a, b))] + offset
            # self.graph[cross_edge[0]][cross_edge[1]]["weight"] = math.inf

    def close_edge(self, u, v):
        '''
        Closes an edge by setting its value to math.inf
        '''
        # print(f"closing edge ({u}-{v})")
        self.graph[u][v]['weight'] = math.inf

    def add_star(self, node_set, weights, pos=None):
        '''
        Adds a single node connected to a given set of vertices (as start for Dijkstra) using the original positions of a node as coordinates
        '''
        a = self.graph.number_of_nodes()
        while a in self.graph:
            a += 1
        if pos != None:
            self.graph.add_node(a, x=pos[0], y=pos[1])
        else:
            self.graph.add_node(a)

        for node in node_set:
            self.graph.add_edge(a, node, weight=weights[node])
        return a

    def remove_star(self, node):
        '''
        Removes a star from the grid (used for stars, which have been added)
        '''
        self.graph.remove_node(node)

    def get_parent(self, node):

        if self.is_port(node):
            return self.port_prnt[node]
        else:
            print(f"{node} is not a port. No parent is defined.")
            input()
            pass