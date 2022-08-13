import math
import os
import sys
import logging

import networkx as nx

import timeit

from ac_metro.exporting.pickleExport import PickleExport
from ac_metro.importing.graphmlImport import GraphmlImport
from ac_metro.render.simpleIpeRender import SimpleIpeRender
from ac_metro.schematic.octiMetroHeuristicSchematic import OctiMetroHeuristicSchematic
from ac_metro.utility.color.colorDictUtility import ColorDictUtility
from ac_metro.utility.custom.octi.helper import bprint, render_grid, render_map, \
    mark_not_valid_candidates, set_meta_data3, is_smooth, is_turning, is_smooth_station
from ac_metro.utility.custom.octi.shapeImport import ShapeImport
from ac_metro.utility.deg2Heuristic.deg2HeuristicUtility import Deg2HeuristicUtility
from ac_metro.utility.deg2Heuristic.undeg2HeuristicUtility import UnDeg2HeuristicUtility
from ac_metro.utility.geometry.averageDistanceUtility import AverageDistanceUtility
from ac_metro.utility.geometry.boundingBoxUtility import BoundingBoxUtility
from ac_metro.utility.geometry.simplificationUtility import Simplifier
from ac_metro.utility.grid.gridCreatorUtility import GridCreator
from ac_metro.utility.grid.gridOverlayCombineUtility import GridOverlayer
from ac_metro.utility.grid.gridRasterizeUtility import GridRasterizer
from ac_metro.utility.paths.spaceAroundPathsUtility import SpaceAroundPathsUtility
from ac_metro.utility.planarization.planarizationUtility import PlanarizationUtility
from ac_metro.utility.planarization.unplanarizationUtility import UnPlanarizationUtility


def mark_smooth_inc_edges(map):
    inc_edges = set()
    for u, v, data in map.graph.edges(data=True):
        if is_smooth(data["edge"]):
            # print(f"for smooth edge {(u, v)}")
            for n in map.graph.neighbors(u):
                # print(f"\tedge {(u, n)} is ", end="")
                if n == v:
                    # print("being looked at")
                    # This is the smooth edge itself
                    continue
                if is_smooth(map.graph[u][n]["edge"]):
                    # print("smooth")
                    # This neighbor is smooth itself
                    continue
                # Now (u, n) is non-smooth and incident
                # print("incident to a smooth edge")
                inc_edges.add((u, n))
                inc_edges.add((n, u))
            for n in map.graph.neighbors(v):
                # print(f"\tedge {(v, n)} is ", end="")
                if n == u:
                    # print("being looked at")
                    # This is the smooth edge itself
                    continue
                if is_smooth(map.graph[v][n]["edge"]):
                    # print("smooth")
                    # This neighbor is smooth itself
                    continue
                # Now (v, n) is non-smooth and incident
                # print("incident to a smooth edge")
                inc_edges.add((v, n))
                inc_edges.add((n, v))

    for u, v in map.graph.edges():
        if (u, v) in inc_edges:
            map.graph[u][v]["smooth_inc"] = True
        else:
            map.graph[u][v]["smooth_inc"] = False


logging.getLogger('matplotlib').setLevel(logging.WARNING)

# Example execution
# python -m src.Tests.octi_metro_heuristic.py M m 1.5 3 20 1 1 1

# Random seed
SEED = 20

# Render Flags
RENDER_INPUT = False
RENDER_BASE_GRID = False
RENDER_BASE_SHAPE = False
RENDER_SIMPLIFIED_SHAPE = False
RENDER_ON_FAIL = False
RENDER_FINAL_GRID = False
RENDER_INPUT_ON_FINAL_GRID = False

# Program
# DEBUG FLAGS
# DEG2 = False                  # SET VIA COMMANDLINE
# PLANARIZE = False             # SET VIA COMMANDLINE
# SPACE_AROUND = False          # SET VIA COMMANDLINE
SHAPE = True

# GRID SETTINGS
# GRID_FACTOR = 0.75            # SET VIA COMMANDLINE
# CAND_SET_FACTOR = 3           # SET VIA COMMANDLINE

# TODO Crossreference with the old paper to validate weights (And check in ExtendedGraph for bend costs)
COST_GRID_EDGE = 20  # TODO: SET VIA COMMANDLINE
COST_GRID_DIAG_EDGE = 20
COST_SHAPE_EDGE = 1  # TODO: SET VIA COMMANDLINE
COST_SHAPE_CONNECTION_EDGE = 10  # TODO: SET VIA COMMANDLINE
COST_TURN_FACTOR = 500  # TODO: SET VIA COMMANDLINE

# SHAPE SETTINGS
UP_DIST_T = 1.5
LOW_DIST_T = 0.5
PROX_REMOVAL = 0.2
# INS_TECH = "distance"
# INS_TECH = "crossings"
SUBSAMPLE_SHAPE = True
SUBSAMPLE_SHAPE_FACTOR = math.sqrt(2)

# ALGORITHM SETTINGS
# MAX_ATTEMPTS = 20             # SET VIA COMMANDLINE
E_ORDER = "random"
# E_ORDER = "line_degree"
# E_ORDER = "length"
PREFER_EDGE_ORDER_ATTRIBUTE = ""
# SKIP_ROUTING = "smooth_inc"
SKIP_ROUTING = ""
ACCEPT_FAILS = 0

# MAP SETTINGS
KEEP_TURNING_STATIONS = True
KEEP_SMOOTH_STATIONS = True
TURN_ATTRIBUTE = "turningpoint"
SMOOTH_ATTRIBUTE = "smoothNode"

# CANDIDATE SET SETTINGS
SKIP_CLOSE_NEIGHBORS = True
SKIP_CLOSE_NEIGHBORS_FACTOR = 0.5
PREFER_CANDIDATE_ATTRIBUTE = "shape_part"  # TODO: move to commandline argument
PREFER_CANDIDATE_FACTOR = 10  # TODO: move to commandline argument
PREFER_CANDIDATE_COUNT = 2

# EXPERIMENTAL FLAGS
SIMPLIFY_SHAPE = False  # unclear
REMOVE_ORIGINAL_SHAPE = False  # unclear
ADD_CROSSING_OF_SHAPE = False  # unclear
LOCAL_SEARCH = int(sys.argv[8]) > 0 if len(sys.argv) > 8 else False  # does not work yet


global_start = timeit.default_timer()

NAME, DIR = set_meta_data3(sys.argv)
print(NAME)
print(DIR)

GRID_FACTOR = float(sys.argv[2])
CAND_SET_FACTOR = float(sys.argv[3])
MAX_ATTEMPTS = int(sys.argv[4])
DEG2 = int(sys.argv[5]) > 0
PLANARIZE = int(sys.argv[6]) > 0
SPACE_AROUND = int(sys.argv[7]) > 0

# PATH_BASE = "/home/soeren/Projects/ac-metro-schematization-framework/ac_metro/tests/Datasets/for_octiSchematic/newnewData"
PATH_BASE = "../../../output-optimisation/"
RUNNING_DIR = PATH_BASE + "/" + DIR + "/"
# OUT_DIR = PATH_BASE + "/OutputFinal/" + f"{DIR}{'-Deg2' if DEG2 else ''}{'-ss'+str(round(SUBSAMPLE_SHAPE_FACTOR, 2)) if SUBSAMPLE_SHAPE else ''}-g{GRID_FACTOR}-c{CAND_SET_FACTOR}-sc{PREFER_CANDIDATE_COUNT}{'-LS' if LOCAL_SEARCH else ''}"
OUT_DIR = PATH_BASE + "/../output-octilinear/" + f"{DIR}-g{GRID_FACTOR}-c{CAND_SET_FACTOR}-p{PROX_REMOVAL}-l{LOW_DIST_T}-u{UP_DIST_T}"

bprint(f"Output is found in {OUT_DIR}")

if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)

# TODO: implement that DEG2 simply retains edges based on contextual information
if DEG2 and PREFER_EDGE_ORDER_ATTRIBUTE != "":
    print("Combination of preferred attribute (e.g. 'smooth') and Deg 2 heuristic does not work yet")

# Reading data
print(f"Reading in the map '{NAME}' from '{RUNNING_DIR}'")
options = {
    "filepath": f"{RUNNING_DIR}output_metro.graphml",
    "x_attr": "x_coordinate",
    "y_attr": "y_coordinate",
    "label_attr": "station_name",
    'line_attr': 'line'
}
if PREFER_EDGE_ORDER_ATTRIBUTE != "":
    options['contextual_edge_information'] = {
        PREFER_EDGE_ORDER_ATTRIBUTE: PREFER_EDGE_ORDER_ATTRIBUTE
    }
options['contextual_station_information'] = {
        TURN_ATTRIBUTE: TURN_ATTRIBUTE,
        SMOOTH_ATTRIBUTE: SMOOTH_ATTRIBUTE
}
metro_map = GraphmlImport.load(options)

# for line, entry in metro_map.getLines():
#     print(line)
#     print(entry)
#     print()

# SimpleIpeRender.render(metro_map, {"filepath": f"{OUT_DIR}", "filename": f"{NAME}_input", "preset": f"{NAME}", "use_labels": False})

renaming_map = {}
for station in metro_map.getStations():
    renaming_map[station.id] = station.label

ColorDictUtility.execute(metro_map, {"preset": NAME})

if RENDER_INPUT:
    print(f"Rendering input")
    render_map(metro_map, nx.Graph(), f"{OUT_DIR}/-1_{DIR}_input.png")

if PLANARIZE:
    print(f"Planarizing input")
    PlanarizationUtility.execute(metro_map)

if DEG2:
    print(f"Using DEG-2 compression")
    deg2HeuristicUtilityOptions = {
        "number_of_nodes_per_path": 0,  # TODO: THIS SHOULD BE 0
        "number_of_nodes_per_path_to_terminal_stations": 0,  # TODO: THIS SHOULD BE 0
        "do_not_add_new_stations": True
    }
    stations_to_keep = []
    for station in metro_map.getStations():
        if KEEP_TURNING_STATIONS and is_turning(station):
            stations_to_keep.append(station.id)
        if KEEP_SMOOTH_STATIONS and is_smooth_station(station):
            stations_to_keep.append(station.id)
    deg2HeuristicUtilityOptions["stations_that_have_to_remain"] = stations_to_keep

    Deg2HeuristicUtility.execute(metro_map, deg2HeuristicUtilityOptions)

    # add smooth information to newly created edges
    for u, v, data in metro_map.graph.edges(data=True):
        smooth_attribute_found = False
        if data["edge"].hasContextualInformation(Deg2HeuristicUtility.UTIL_NAME()):
            removedEdges = data["edge"].getContextualInformation(Deg2HeuristicUtility.UTIL_NAME()).removedEdges
            for removedEdge in removedEdges:
                removedEdgeData = removedEdge['edge']
                if is_smooth(removedEdgeData):
                    smooth_attribute_found = True
                    break
            data["edge"].setContextualInformation(GraphmlImport.IMPORT_NAME(), {"smooth": smooth_attribute_found})

# Marking edge, which are incident to smooth edges
print(f"Marking non-smooth edges incident to smooth edges")
mark_smooth_inc_edges(metro_map)

# Preparing the grid
print("Rastering the map")
# Computing the dimensions of the input map
avg_dist = AverageDistanceUtility.execute(metro_map, {})

# Compute cell size of the grid based on map dimensions and GRID_FACTOR
grid_cell_size = GRID_FACTOR * avg_dist

# compute number of grid columns and rows, based on map dimensions
bb_min_x, bb_min_y, bb_width, bb_height = BoundingBoxUtility.execute(metro_map, {})
grid_width, grid_height = GridRasterizer.execute(metro_map, {"resolution": grid_cell_size, "bb_width": bb_width, "bb_height": bb_height})


print("Creating grid")
# Creating the grid
grid_options = {
    "width": grid_width,
    "height": grid_height,
    "x_scale": grid_cell_size,
    "y_scale": grid_cell_size,
    "bottom_left": (bb_min_x, bb_min_y),
    "type": "octolinear",
    "grid_cost": COST_GRID_EDGE,
    "grid_diagonal_cost": COST_GRID_DIAG_EDGE
}
grid_graph = GridCreator.execute(metro_map, grid_options)

if RENDER_INPUT:
    print(f"Rendering input on grid")
    render_map(metro_map, grid_graph, f"{OUT_DIR}/0_{DIR}_input_on_grid.png")

if RENDER_BASE_GRID:
    print(f"Rendering base grid")
    render_grid(grid_graph, f"{OUT_DIR}/1_{DIR}_base_grid.png")

if SHAPE:
    # Reading in the shape, placing it on the map and intersects it into the grid (simplifyication of the shape can
    # be done)
    print(f"Reading in the shape, ID starts at {grid_graph.number_of_nodes()}")
    options = {
        "filepath": f"{RUNNING_DIR}/output_guide.graphml",  # TODO should be variable based on commandline argument
        "x_attr": "x_coordinate",
        "y_attr": "y_coordinate",
        "id_start": grid_graph.number_of_nodes(),
        "shape_cost": COST_SHAPE_EDGE
    }
    shape = ShapeImport.load(options)

    options = {
        "filepath": f"{RUNNING_DIR}/output_guide.graphml",  # TODO should be variable based on commandline argument
        "x_attr": "x_coordinate",
        "y_attr": "y_coordinate",
        "label_attr": "station_name",
        'line_attr': 'line'
    }
    metro_map_2 = GraphmlImport.load(options)
    # SimpleIpeRender.render(metro_map_2, {"filepath": f"{OUT_DIR}", "filename": f"{NAME}_shape", "use_labels": False})


    if RENDER_BASE_SHAPE:
        print(f"Rendering shape")
        render_grid(shape, f"{OUT_DIR}/3_{NAME}_just_shape.png")

    if SUBSAMPLE_SHAPE:
        id_start = grid_graph.number_of_nodes() + shape.number_of_nodes()
        new_id = id_start
        print(f"Subsampling shape")
        nodes_to_add = {}
        edge_to_add = []
        edge_to_remove = []
        for u, v in shape.edges():
            # print(f"\tconsidering {(u, v)}")
            xu = shape.nodes()[u]['x']
            yu = shape.nodes()[u]['y']
            xv = shape.nodes()[v]['x']
            yv = shape.nodes()[v]['y']
            distance = math.dist((xu, yu), (xv, yv))
            # print(f"\tdistance: {distance} and threshhold is {(grid_cell_size * SUBSAMPLE_SHAPE_FACTOR)} result: {distance/(grid_cell_size * SUBSAMPLE_SHAPE_FACTOR)}")
            number = math.floor(distance/(grid_cell_size * SUBSAMPLE_SHAPE_FACTOR))
            if number >= 1:
                # print(f"\tAdding {number} vertices to {(u, v)}")
                edge_to_remove.append((u, v))
                last = u
                for i in range(1, number+1):
                    ratio = float(i)/float(number+1)
                    edge_to_add.append((last, new_id))
                    nodes_to_add[new_id] = (xu * (1-ratio) + xv * (ratio), yu * (1-ratio) + yv * (ratio))
                    last = new_id
                    new_id += 1
                edge_to_add.append((new_id-1, v))

        for node in nodes_to_add:
            coords = nodes_to_add[node]
            shape.add_node(node, x=coords[0], y=coords[1])
        for edge in edge_to_add:
            shape.add_edge(edge[0], edge[1], crossing=[], weight=COST_SHAPE_EDGE)
            # print(f"added {edge}")
        for edge in edge_to_remove:
            shape.remove_edge(edge[0], edge[1])

    if SIMPLIFY_SHAPE:
        print(f"Simplifying the shape")
        shape = Simplifier.execute(None, {"graph": shape, "factor": grid_cell_size / 10,
                                          "id_start": grid_graph.number_of_nodes()})

    if RENDER_BASE_SHAPE:
        print(f"Rendering shape")
        render_grid(shape, f"{OUT_DIR}/4_{NAME}_final_shape.png")

    print(f"Combining grid and shape using distance method")
    overlay_options = {
        "grid": grid_graph,
        "shape": shape,
        "upper_threshhold": UP_DIST_T * grid_cell_size,
        "lower_threshhold": LOW_DIST_T * grid_cell_size,
        "shape_connection_cost": COST_SHAPE_CONNECTION_EDGE,
        "proximity_removal": PROX_REMOVAL * grid_cell_size,
        "method": "distance"
    }

    custom_grid = GridOverlayer.execute(metro_map, overlay_options)

    # make sure too close nodes are not all valid candidate nodes
    if SKIP_CLOSE_NEIGHBORS:
        print(f"Skipping close neighbors")
        mark_not_valid_candidates(custom_grid, grid_cell_size * SKIP_CLOSE_NEIGHBORS_FACTOR)

else:
    print("No shape added. Grid is used as is")

    # We mark nodes not on the shape as such
    for node, data in grid_graph.nodes(data=True):
        grid_graph.nodes()[node]["shape_part"] = False

    custom_grid = grid_graph

if RENDER_FINAL_GRID:
    print(f"Rendering final grid")
    render_grid(custom_grid, f"{OUT_DIR}/5_{DIR}_final_grid.png")

if RENDER_INPUT_ON_FINAL_GRID:
    print(f"Rendering input")
    render_map(metro_map, custom_grid, f"{OUT_DIR}/6_{DIR}_input_final.png")

print("Starting schematization")
# Schematizing
options = {
    "name": NAME,
    "run_dir": RUNNING_DIR,
    "out_path_dir": OUT_DIR,
    "method": "PortNodes",
    "max_attempts": MAX_ATTEMPTS,
    "grid": custom_grid,
    "cell_size": grid_cell_size,
    "bend_cost_factor": COST_TURN_FACTOR,  # DONE: check how this is handled internally in ExtendedeGrid.
    "edge_order_method": E_ORDER,
    "cand_set_method": "circle",
    "cand_set_factor": CAND_SET_FACTOR,
    "candidate_ignore": "not_candidate",
    "break_on_failure": True,
    "local_search_optimization": LOCAL_SEARCH,
    "render": False,
    "crossing_variant": 3
}
if SEED > 0:
    options["seed"]: SEED
if ACCEPT_FAILS > 0:
    options["accept_failed_edges"] = ACCEPT_FAILS
if SKIP_ROUTING != "":
    options["skip_routing"] = SKIP_ROUTING
if PREFER_EDGE_ORDER_ATTRIBUTE != "":
    options["prefer_in_edge_order"] = PREFER_EDGE_ORDER_ATTRIBUTE
if PREFER_CANDIDATE_ATTRIBUTE != "":
    options["prefer_as_candidate"] = PREFER_CANDIDATE_ATTRIBUTE
    options["prefer_as_candidate_factor"] = PREFER_CANDIDATE_FACTOR
    options["prefer_cand_count"] = PREFER_CANDIDATE_COUNT


grid, success = OctiMetroHeuristicSchematic.schematize(metro_map, options)

# for u, v, data in metro_map.graph.edges(data=True):
#     print(f"{(u, v)} has: ", end="")
#     print(data["edge"].getControlPoints())


# Drawing the result
success_name = "success" if success else "failure"
if DEG2:
    print("Reverting DEG-2 compression")
    UnDeg2HeuristicUtility.execute(metro_map, {"preserve_control_points": True})
if PLANARIZE:
    print("Reverting planarization")
    UnPlanarizationUtility.execute(metro_map)

# for ID in renaming_map:
#     print(ID)
#     renaming_station = metro_map.getStation(ID)
#     renaming_station.label = renaming_map[ID]


# render_map(metro_map, grid, f"{OUT_DIR}/9_{DIR}_final.png")



# PickleExport.export(metro_map, {"filepath": f"{OUT_DIR}/9_{NAME}_final.map"})

if SPACE_AROUND:
    print("Spacing out stations on paths")
    SpaceAroundPathsUtility.execute(metro_map)
    render_map(metro_map, grid_graph, f"{OUT_DIR}/10_{DIR}_after_spacing.png")

global_end = timeit.default_timer()

print('Global Time: ', global_end - global_start)

# PickleExport.export(metro_map, {"filepath": f"{OUT_DIR}/10_{DIR}_after_spacing.map"})
SimpleIpeRender.render(metro_map, {"filepath": f"{OUT_DIR}", "filename": f"{NAME}", "use_labels": True})
