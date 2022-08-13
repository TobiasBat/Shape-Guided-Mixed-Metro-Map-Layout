from shapely.speedups._speedups import LineString

from ac_metro.utility.abstractUtility import AbstractUtility

import bentley_ottmann.planar as bent_ott
import ground.base as g_base
from ac_metro.utility.custom.octi.poly_point_isect import *


class GridCrossing(AbstractUtility):

    @staticmethod
    def execute(map, options=None):
        grid = options["grid"]
        if "variant" in options and options["variant"] == 2:
            return compute_crossings_in_grid_2(grid)
        elif "variant" in options and options["variant"] == 3:
            return compute_crossings_shape_grid(grid)
        else:
            return compute_crossings_in_grid(grid)


def compute_crossings_in_grid(grid):
    edge_segments, seg_to_edge = convert_edges_to_segments(grid)
    intersections = bent_ott.segments_intersections(edge_segments)
    for int in intersections:
        edge1, edge2 = seg_to_edge[int[0]], seg_to_edge[int[1]]
        u1, v1, u2, v2 = edge1[0], edge1[1], edge2[0], edge2[1]
        if u1 == u2 or u1 == v2 or v1 == u2 or v1 == v2:
            continue
        grid[u1][v1]["crossing"].append(edge2)
        grid[u2][v2]["crossing"].append(edge1)
    return grid


def convert_edges_to_segments(grid):
    # This is following an example from the pypi website. Not entirely sure whats happening with the context
    context = g_base.get_context()
    Point, Segment = context.point_cls, context.segment_cls
    segments = []
    segment_id_to_edge = {}
    ID = 0
    for edge in grid.edges():
        id1, id2 = edge[0], edge[1]
        node1, node2 = grid.nodes[id1], grid.nodes[id2]
        # print(f"{node1}\n{node2}\n")
        segment = Segment(Point(node1['x'], node1['y']), Point(node2['x'], node2['y']))
        segments.append(segment)
        segment_id_to_edge[ID] = edge
        ID += 1
    return segments, segment_id_to_edge


def compute_crossings_in_grid_2(grid):
    segments = []
    coords_to_edge = {}
    for edge in grid.edges():
        id1, id2 = edge[0], edge[1]
        node1, node2 = grid.nodes[id1], grid.nodes[id2]
        segment = ((node1['x'], node1['y']), (node2['x'], node2['y']))
        segments.append(segment)
        coords_to_edge[((node1['x'], node1['y']), (node2['x'], node2['y']))] = edge
        coords_to_edge[((node2['x'], node2['y']), (node1['x'], node1['y']))] = edge
    crossings = isect_segments_include_segments(segments)
    for int in crossings:
        seg1, seg2 = int[1][0], int[1][1]
        edge1, edge2 = coords_to_edge[seg1], coords_to_edge[seg2]
        u1, v1, u2, v2 = edge1[0], edge1[1], edge2[0], edge2[1]
        if u1 == u2 or u1 == v2 or v1 == u2 or v1 == v2:
            print("should never happen!")
        grid[u1][v1]["crossing"].append(edge2)
        grid[u2][v2]["crossing"].append(edge1)
    return grid

def compute_crossings_shape_grid(grid):

    for edge in grid.edges():
        if edge[0] == -1 or edge[1] == -1:
            print(edge)
            input()
        if not grid[edge[0]][edge[1]]["shape_part"]:
            continue
        print(grid.nodes()[edge[0]])
        x1 = grid.nodes()[edge[0]]['x']
        y1 = grid.nodes()[edge[0]]['y']
        x2 = grid.nodes()[edge[1]]['x']
        y2 = grid.nodes()[edge[1]]['y']
        line = LineString([(x1, y1), (x2, y2)])
        for edge2 in grid.edges():
            if edge2[0] == -1 or edge2[1] == -1:
                print(edge2)
                input()
            if grid[edge2[0]][edge2[1]]["shape_part"]:
                continue
            x3 = grid.nodes()[edge2[0]]['x']
            y3 = grid.nodes()[edge2[0]]['y']
            x4 = grid.nodes()[edge2[1]]['x']
            y4 = grid.nodes()[edge2[1]]['y']
            line2 = LineString([(x3, y3), (x4, y4)])
            if line.intersects(line2):
                grid[edge[0]][edge[1]]["crossing"].append(edge2)
                grid[edge2[0]][edge2[1]]["crossing"].append(edge)
    return grid
