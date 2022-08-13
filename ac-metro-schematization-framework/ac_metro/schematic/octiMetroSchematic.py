# NEVER FINISHED CLASS!
import itertools
import math

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from ac_metro.map.map import Map
from ac_metro.schematic.abstractSchematic import AbstractSchematic

hard_edge_costs = {
    "sink": 1 + math.pi,
    "jump": 1
}
bend_cost_offset = 0

class GridNode:
    def __init__(self, id, x=None, y=None, is_port=False, parent_id=None):
        self.id = id
        self.x = x
        self.y = y
        self.is_port = is_port
        self.inc_edges = []
        self.angles = {}
        if is_port:
            self.parent = parent_id

    def addEdge(self, edge):
        '''
        This sorts the added edges based on the atan2 centered at this station
        '''
        other = edge.source
        if other == self:
            other = edge.target
        angle = self.angleTo(other)
        self.angles[other.id] = angle
        for i in range(len(self.inc_edges)):
            this_edge = self.inc_edges[i]
            this_other = this_edge.source
            if this_other == self:
                this_other = this_edge.target
            if angle < angle[this_other.id]:
                self.inc_edges.insert(i, edge)

    def angleTo(self, other):
        return math.atan2(other.y - self.y, other.x - self.x)


class GridEdge:
    def __init__(self, source, target, type, cost):
        # Types are:
        # 1: port to parent
        # 2: port to port (intra node)
        # 3: port to port (inter node)
        self.source = source
        self.target = target
        self.type = type
        self.cost = cost
        self.crossing = set()

    def addCrossingEdge(self, edge):
        self.crossing.add(edge)

    def crosses(self, edge):
        return edge in self.crossing


def on_segment(p, q, r):
    if r[0] <= max(p[0], q[0]) and r[0] >= min(p[0], q[0]) and r[1] <= max(p[1], q[1]) and r[1] >= min(p[1], q[1]):
        return True
    return False


def orientation(p, q, r):
    val = ((q[1] - p[1]) * (r[0] - q[0])) - ((q[0] - p[0]) * (r[1] - q[1]))
    if val == 0 : return 0
    return 1 if val > 0 else -1


class Grid:
    def __init__(self):
        self.grid = nx.Graph()

    def get_neighbors(self, node_id):
        self.grid

    def addNode(self, g_node: GridNode):
        self.grid.add_node(g_node.id, g_node=g_node)

    def get_nodes(self):
        nodeSet = []
        for node_id, data in self.grid.nodes(data=True):
            nodeSet.append(data['g_node'])
        return nodeSet

    def get_pos_dict(self):
        pos = {}
        for node in self.get_nodes():
            pos[node.id] = np.array([node.x, node.y])
        return pos

    def addConnection(self, g_edge: GridEdge):
        for source, target, data in self.grid.edges(data=True):
            existing_edge = data["g_edge"]
            if existing_edge.type != 3:
                continue
            if self.edges_cross(existing_edge, g_edge):
                g_edge.addCrossingEdge(existing_edge)
                existing_edge.addCrossingEdge(g_edge)
        self.grid.add_edge(g_edge.source, g_edge.target, g_edge=g_edge)

    def getAngleBetween(self, id1, id2, id3):
        n1, n2, n3 = self.grid.nodes[id1]["g_node"], self.grid.nodes[id2]["g_node"], self.grid.nodes[id3]["g_node"]
        ang = math.atan2(n3.y - n2.y, n3.x - n2.x) - math.atan2(n2.y - n1.y, n2.x - n1.x)
        return ang + math.pi if ang < 0 else ang

    def edges_cross(self, edge1, edge2):
        s1, t1, s2, t2 = edge1.source, edge1.target, edge2.source, edge2.target
        p1 = (self.grid.nodes[s1]["g_node"].x, self.grid.nodes[s1]["g_node"].y)
        q1 = (self.grid.nodes[t1]["g_node"].x, self.grid.nodes[t1]["g_node"].y)
        p2 = (self.grid.nodes[s2]["g_node"].x, self.grid.nodes[s2]["g_node"].y)
        q2 = (self.grid.nodes[t2]["g_node"].x, self.grid.nodes[t2]["g_node"].y)

        o1 = orientation(p1, q1, p2)
        o2 = orientation(p1, q1, q2)
        o3 = orientation(p2, q2, p1)
        o4 = orientation(p2, q2, q1)

        if o1 != o2 and o3 != o4:
            return True

        if o1 == 0 and on_segment(p1, q1, p2): return True
        if o2 == 0 and on_segment(p1, q1, q2): return True
        if o3 == 0 and on_segment(p2, q2, p1): return True
        if o4 == 0 and on_segment(p2, q2, q1): return True
        return False


def angleBendCost(grid, node_id, port_ID_1, port_ID_2):
    return 1 + bend_cost_offset + (math.pi - grid.getAngleBetween(port_ID_1, node_id, port_ID_2))


def addDetailToGrid(node_set, edge_set, calculateBendCost):
    '''
    Given a set of nodes and edges and a function to determine the cost of edges, this creates the necessary detail grid
    '''
    node_count = len(node_set)
    grid = Grid()

    # incidence dict gives us info about the connections and edges
    incidences = {}
    for edge_tuple in edge_set:
        a, b = edge_tuple[0], edge_tuple[1]
        if a not in incidences:
            incidences[a] = []
        if b not in incidences:
            incidences[b] = []
        incidences[a].append(b)
        incidences[b].append(a)

    ID = node_count
    port_connections = {} #maps tuple (a,b) for edge (a,b) to a tuple of port ids (c,d), so we know which ports to connect in the inter connection step
    for node_id, coords in node_set.items():
        x, y = coords[0], coords[1]
        start_ID = ID
        grid.addNode(GridNode(node_id, x, y, is_port=False, parent_id=None))
        incident_edges = incidences[node_id]
        for target_ID in incident_edges:
            grid.addNode(GridNode(ID, x, y, is_port=True, parent_id=node_id))
            if frozenset((node_id, target_ID)) not in port_connections:
                port_connections[frozenset((node_id, target_ID))] = [] #Next time we see it the roles will be reversed
            port_connections[frozenset((node_id, target_ID))].append(ID) #(node_id, target_id) in any order now maps to a list with two ids, which are the corresponding ports
            ID += 1
        for port_ID in range(start_ID, ID):
            grid.addConnection(GridEdge(node_id, port_ID, 1, hard_edge_costs["sink"]))
        for port_ID_1, port_ID_2 in itertools.combinations(range(start_ID, ID), 2):
            grid.addConnection(GridEdge(port_ID_1, port_ID_2, 2, calculateBendCost(grid, node_id, port_ID_1, port_ID_2)))
        # Metanode connections are set up
    # All Metanodes are added (including internal edges)

    for edge_tuple in edge_set:
        a, b = edge_tuple[0], edge_tuple[1]
        c, d = port_connections[frozenset((a, b))][0], port_connections[frozenset((a, b))][1]
        edge_to_add = GridEdge(c, d, 3, hard_edge_costs["jump"])
        grid.addConnection(edge_to_add)

    return grid


def createRegularGrid(width, height, x_scale, y_scale):
    ID = 0
    node_set, edge_set = {}, set()
    for i in range(height):
        for j in range(width):
            node_set[ID] = (i*x_scale, j*y_scale)
            if j > 0:
                edge_set.add((ID, ID-1))
            if i > 0:
                edge_set.add((ID, ID-width))
            if i > 0 and j > 0:
                edge_set.add((ID, ID-width-1))
            if i > 0 and j < width:
                edge_set.add((ID, ID-width+1))
            ID += 1
    return addDetailToGrid(node_set, edge_set, angleBendCost)

class OctiMetroSchematic(AbstractSchematic):

    @staticmethod
    def schematize(map: Map, options=None):
        grid = createRegularGrid(3, 3, 10, 10)
        print(pos)
        nx.draw(grid.grid, pos, with_labels = True)#, cmap=plt.get_cmap('jet'), node_color=values, node_size=500)
        plt.show()
        print(f"done")
