import timeit
from cProfile import label
import copy
import math

import json

import networkx as nx
import numpy as np
from ac_metro.exporting.pickleExport import PickleExport
from ac_metro.importing.graphmlImport import GraphmlImport
from ac_metro.map.map import Map
from ac_metro.render.matplotlibGridRender import MatplotlibGridRender
from ac_metro.render.simpleIpeRender import SimpleIpeRender
from ac_metro.schematic.abstractSchematic import AbstractSchematic
from ac_metro.utility.custom.octi.extendedGrid import ExtendedGrid
from ac_metro.utility.custom.octi.helper import dist, is_smooth, render_map, bprint, is_smooth_station, render_grid
from ac_metro.utility.geometry.circularVertexOrderUtility import CircularVertexOrderUtility
from ac_metro.utility.custom.octi.edgeOrderUtility import EdgeOrderUtility
from ac_metro.utility.grid.gridCrossingsUtility import GridCrossing
from ac_metro.utility.grid.gridPlacementUtility import GridPlacementUtility
from ac_metro.utility.grid.mapCenteringUtility import MapCenterer


class OctiMetroHeuristicSchematic(AbstractSchematic):
    NAME = None
    OUT_DIR = None
    settled_edges_costs = None
    settled_edges = None
    settled_nodes = None
    orig_positions = None

    @staticmethod
    def schematize(map: Map, options=None):
        OctiMetroHeuristicSchematic.RUNNING_DIR = options["run_dir"] if "run_dir" in options else "NODIR"
        OctiMetroHeuristicSchematic.NAME = options["name"] if "name" in options else "NONAME"
        OctiMetroHeuristicSchematic.OUT_DIR = options["out_path_dir"] if "out_path_dir" in options else "DEFAULT"

        if "grid" in options:
            grid = options["grid"]
        else:
            print("Please supply a base-grid.")
            return

        print("Adding grid properties")
        for u, v in grid.edges():
            grid[u][v]["crossing"] = []

        print("Computing crossing in grid")
        grid = GridCrossing.execute(None, {"grid": grid})

        print("Centering map on grid")
        # Centering map on the grid
        # MapCenterer.execute(map, {"grid": grid})

        options["grid"] = grid

        return OctiMetroHeuristicSchematic.schematize_port_nodes(map, options)

    @classmethod
    def schematize_port_nodes(cls, map: Map, options=None):

        print(options)

        MAX_ATTEMPTS = options["max_attempts"] if "max_attempts" in options else 1
        RENDER = "render" in options and options["render"]
        simple_grid = options["grid"]
        # grid_cell_size = options["cell_size"]
        SKIP_ROUTING = ""
        if "skip_routing" in options:
            SKIP_ROUTING = options["skip_routing"]

        print(f"Computing extended grid")
        grid = ExtendedGrid(simple_grid, options['bend_cost_factor'])

        print(f"Extended grid has {max(grid.graph.nodes())}")

        ###################### CODE ADDED TO UNDERSTAND THE EXTENDED GRID #######################
        # also numpy import 

        H = nx.subgraph(grid.graph, [0,1982,1984,1980,1,1981,1986,1988,1990,1992,45,2336, 2338, 2334, 1989, 1983, 46, 1985, 2335, 2342, 2344, 2346, 2340, 1997, 1991])
         
        pos = {}
        for node, data in H.nodes(data=True):
            # print(f"DEBUG: {node} has {data}")
            pos[node] = np.array([data["x"], data["y"]])

        labels = nx.get_edge_attributes(H,'weight')
        for label,value in labels.items():
            labels[label]=np.round(value,3)

        #nx.draw(H,pos, with_labels=True, font_weight='bold', node_color='lightblue', node_size=200,connectionstyle='arc3, rad = 0.1')
        #nx.draw_networkx_edge_labels(H,pos,edge_labels=labels)


        # if RENDER:
        # print(f"Rendering extended grid")
        # render_grid(grid.graph, f"{cls.OUT_DIR}/{cls.NAME}_extended_grid.png")

        print("Saving original positions")
        # We keep the original positions of all stations in the map
        cls.orig_positions = {}
        for node, data in map.graph.nodes(data=True):
            cls.orig_positions[node] = (data["station"].x, data["station"].y)

        # Sanity check for shape_marking
        if "prefer_as_candidate" in options:
            prefer_cand_att = options["prefer_as_candidate"]
            for node, data in grid.graph.nodes(data=True):
                if prefer_cand_att not in data and not grid.is_port(node):
                    print(f"{prefer_cand_att} should be preferred as candidate, but not all grid nodes have this tag")

        # We compute the circular order around every station in the map
        # TODO: Check this for correctness
        print(f"Computing circular vertex order")
        map_circ_ord = CircularVertexOrderUtility.execute(map, {"min_degree": 1})

        failure = True
        attempts = 0

        start = timeit.default_timer()

        print(f"Starting routing attempt {attempts + 1}")
        while failure and attempts < MAX_ATTEMPTS:

            if attempts > 0 and "accept_failed_edges" in options and options["accept_failed_edges"] >= len(failed_edges):
                break
            # Resetting to non-routed edges
            map.reset_control_points()
            grid.reset_edge_weights()
            map.set_station_positions(cls.orig_positions)

            # Computing the candidate positions
            print(f"Computing candidate sets")
            cls.candidates = cls.compute_cand_poss(map, options)
            # print(cls.candidates)
            # input()

            # for node in map.graph.nodes():
            #     print(f"cand_pos: {cand_pos[node]}")
            #     print(f"candidates: {candidates[node]}")

            edges_routed = 0  # count of edges attempted to route
            cls.settled_nodes = {}  # ID -> ID in grid
            cls.settled_edges = {}  # (ID1, ID2) -> chain of IDs in grid (polyline)
            cls.settled_edges_costs = {}  # (ID1, ID2) -> int (cost of polyline)
            failed_edges = set()  # failed to route these edges
            failure = False  # At least one failed edge (Redundant?)

            # Compute edge order TODO: figure out if this belongs here or in a Utility

            if attempts > 0:
                print(f"Shuffling edges...")
            edge_order_options = {
                "method": options["edge_order_method"] if "edge_order_method" in options and attempts <= 1 else "random"
            }
            if "seed" in options:
                edge_order_options["seed"] = options["seed"]
            if "prefer_in_edge_order" in options:
                edge_order_options["prefer_att"] = options["prefer_in_edge_order"]

            if "prefer_in_edge_order" in options:
                edge_list_a = []
                edge_list_b = []
                for u, v, data in map.graph.edges(data=True):
                    # print(f"considering{(u, v)}")
                    if data["edge"].hasContextualInformation(GraphmlImport.IMPORT_NAME()) \
                            and "smooth" in data["edge"].getContextualInformation(GraphmlImport.IMPORT_NAME()) \
                            and data["edge"].getContextualInformation(GraphmlImport.IMPORT_NAME())["smooth"]:
                        edge_list_a.append((u, v))
                    else:
                        edge_list_b.append((u, v))

                # print(f"options are {edge_order_options}")
                edge_order_options["edge_list"] = edge_list_a
                # print(f"edgelist is first:\n{edge_order_options['edge_list']}")
                ordered_edges = EdgeOrderUtility.execute(map, edge_order_options)
                # print(f"ordered edges are first:\n{ordered_edges}")

                edge_order_options["edge_list"] = edge_list_b
                # print(f"edgelist is second:\n{edge_order_options['edge_list']}")
                ordered_edges.extend(EdgeOrderUtility.execute(map, edge_order_options))
                # print(f"ordered edges are second:\n{ordered_edges}")

            else:
                edge_order_options["edge_list"] = [(u, v) for u, v in map.graph.edges()]
                ordered_edges = EdgeOrderUtility.execute(map, edge_order_options)

            for id1, id2 in ordered_edges:

                value1, value2 = int(id1[1:]), int(id2[1:])
                id1 = "n" + str(min(value1, value2))
                id2 = "n" + str(max(value1, value2))

                # if id2 == map.getEdge(id1, id2).source:
                #     idtemp = id2
                #     id2 = id1
                #     id1 = idtemp

                # Sanity Check
                if (id1, id2) in failed_edges:
                    print("should not happen")
                    continue
                # TODO: This should rather be something along the lines of an option["dont_route"] = set of edges
                if SKIP_ROUTING != "" and SKIP_ROUTING in map.graph[id1][id2] and map.graph[id1][id2][SKIP_ROUTING]:
                    print(f"skipping {(id1, id2)} because it is {SKIP_ROUTING}")
                    continue

                print(f"\tRouting: ({id1} - {id2})")

                # This is where the magic happens. TODO: CHECK IN DETAIL!!!
                path = cls.route_single_edge(map, grid, simple_grid, id1, id2, map_circ_ord, cls.candidates)

                # Computing the cost of the found path now with the temporary nodes removed
                path_cost = 0
                for i in range(1, len(path)):
                    path_cost += grid.graph[path[i - 1]][path[i]]["weight"]

                # If we did not find a valid path or only one, which uses a closed edge and therefore has cost
                # infinity, we break
                if len(path) == 0 or math.isinf(path_cost):
                    failure = True
                    print(f"\tERROR: Could not route {edges_routed}-th edge from {id1} to {id2}.")
                    failed_edges.add((id1, id2))
                    if "break_on_failure" in options and options["break_on_failure"]:
                        break
                    if "accept_failed_edges" in options and options["accept_failed_edges"] < len(failed_edges):
                        print(options["accept_failed_edges"])
                        print(len(failed_edges))
                        break
                else:
                    # Keep track of settled nodes and edges
                    cls.settled_nodes[id1] = path[0]
                    cls.settled_nodes[id2] = path[-1]

                    cls.set_settled_edge(id1, id2, path)
                    # cls.settled_edges[(id1, id2)] = path
                    # cls.settled_edges[(id2, id1)] = path[::-1]

                    cls.set_settled_edge_cost(id1, id2, path_cost)
                    # cls.settled_edges_costs[(id1, id2)] = path_cost
                    # cls.settled_edges_costs[(id2, id1)] = path_cost

                    # Close edges and nodes, which were used by this path
                    # print(f"Current path: {path}")

                    ########### DRAW TO SEE WHAT IS HAPPENING IN THE EDGES AFTER THE PATH IS SETTLED ###############


                    grid.close_path(path)
                    grid.close_node(path[0])
                    grid.close_node(path[-1])
                    cls.draw_network_state(grid.graph,path)

                    # update station positions
                    station1, station2 = map.graph.nodes[id1]["station"], map.graph.nodes[id2]["station"]
                    station1.x, station1.y = grid.coord(path[0])
                    station2.x, station2.y = grid.coord(path[-1])


                    # if id1 == "n1" or id1 == "n84":
                    #     if id2 == "n1" or id2 == "n84":
                    #         print(f"edge {(id1, id2)} has source {map.getEdge(id1, id2).source} and target {map.getEdge(id1, id2).target}")
                    #         input()


                    # update edge control points
                    ## ASK why here the step is two 
                    control_points = []
                    for node in path[2:-2:2]:
                        # print('control_points :'+ str(node))
                        control_points.append(grid.coord(node))

                    # if id1 == map.getEdge(id1, id2).source:
                    #     for node in path[2:-2:2]:
                    #         control_points.append(grid.coord(node))
                    # else:
                    #     for node in path[-2:2:-2]:
                    #         control_points.append(grid.coord(node))

                    map.getEdge(id1, id2).setControlPoints(control_points, start=id1)
                    # map.getEdge(id1, id2).setControlPoints(control_points, start=id1)

                edges_routed += 1
            attempts += 1

            pickle_options = {
                "filepath": f"{cls.OUT_DIR}/7_{cls.NAME}_{'fail' if failure else 'success'}_{attempts}_{len(failed_edges)}.map"
            }
            # PickleExport.export(map, pickle_options)

            if RENDER:
                print(f"Failed to route {len(failed_edges)} edges")

                SimpleIpeRender.render(map, {"filepath": f"{cls.OUT_DIR}", "filename": f"7_{cls.NAME}_{'fail' if failure else 'success'}_{attempts}_{len(failed_edges)}"})
                render_map(map, grid.graph, f"{cls.OUT_DIR}/7_{cls.NAME}_{'fail' if failure else 'success'}_{attempts}_{len(failed_edges)}.png")
                # MatplotlibGridRender.render(map, {
                #     "grid": grid.graph,
                #     "just_grid": False,
                #     "grid_labels": False,
                #     "show_nodes": False,
                #     "station_names": True,
                #     "draw_straight": SKIP_ROUTING,
                #     "to_file": f"{cls.OUT_DIR}/{cls.NAME}_{attempts}_{'fail' if failure else 'success'}_{len(failed_edges)}.png"
                # })

            # TODO: This NEEDS to be sorted out ASAP
            if not failure and "local_search_optimization" in options and options["local_search_optimization"]:
                #if RENDER:
                #render_map(map, grid.graph, f"{cls.OUT_DIR}/8_{cls.NAME}_before_LS_{attempts}.png")
                cls.perform_local_search_alt(map, grid, simple_grid, cls.candidates, map_circ_ord, attempts)
                #if RENDER:
                #render_map(map, grid.graph, f"{cls.OUT_DIR}/8_{cls.NAME}_after_LS_{attempts}.png")


        stop = timeit.default_timer()
        print('Time: ', stop - start)
        return grid.graph, not failure#, cls, grid, map_circ_ord

    @classmethod
    def route_single_edge(cls, map, grid, simple_grid, id1, id2, map_circ_ord, candidates, opening=True):
        # Variable setup
        station1, station2 = map.graph.nodes[id1]["station"], map.graph.nodes[id2]["station"]

        # If a station is already settled, we know its porition, otherwise we need to consider its candidate set
        src_pos = cls.settled_nodes[id1] if id1 in cls.settled_nodes else -1
        tgt_pos = cls.settled_nodes[id2] if id2 in cls.settled_nodes else -1
        # print(f"srcpos: {src_pos}")
        # print(f"tgtpos: {tgt_pos}")

        added_src, added_tgt = False, False

        # print(f"when routing candidates is {[entry[0] for entry in candidates[id1]]} and cand_pos is {cand_pos[id1]} which is equal {[entry[0] for entry in candidates[id1]] == cand_pos[id1]}")
        # if [entry[0] for entry in candidates[id1]] != cand_pos[id1]:
        #     print("\tBANG! 1")
        # print(f"when routing candidates is {[entry[0] for entry in candidates[id2]]} and cand_pos is {cand_pos[id2]} which is equal {[entry[0] for entry in candidates[id2]] == cand_pos[id2]}")
        # if [entry[0] for entry in candidates[id2]] != cand_pos[id2]:
        #     print("\tBANG! 2")

        cand_pos_src = [entry[0] for entry in candidates[id1]]
        cand_pos_tgt = [entry[0] for entry in candidates[id2]]

        if src_pos == -1 and tgt_pos == -1:
            # print("None settled")
            # removing src, tgt candidate overlap since neither are settled
            # cand_positions1, cand_positions2 = cls.remove_candidate_overlap(cand_pos[id1], cand_pos[id2],
            #                                                                 (station1.x, station1.y),
            #                                                                 (station2.x, station2.y), simple_grid)
            cand_positions1, cand_positions2 = cls.remove_candidate_overlap(cand_pos_src, cand_pos_tgt,
                                                                            (station1.x, station1.y),
                                                                            (station2.x, station2.y), simple_grid)

            # Setting offset weights for candidate source positions
            src_pos = OctiMetroHeuristicSchematic.prepare_unsettled_node(candidates, id1, grid, cand_positions1)
            added_src = True

            # Setting offset weights for candidate target positions
            tgt_pos = OctiMetroHeuristicSchematic.prepare_unsettled_node(candidates, id2, grid, cand_positions2)
            added_tgt = True

        elif tgt_pos == -1:
            # print("target not settled")
            # removing src, from tgt candidates
            # cand_positions2 = cand_pos[id2]
            cand_positions2 = cand_pos_tgt
            if src_pos in cand_positions2:
                # print(f"before removing {src_pos} cand_pos is\n{cand_pos[id2]}\nand candidates is\n{candidates[id2]}")
                cand_positions2.remove(src_pos)
                remove_index = 0
                for entry in candidates[id2]:
                    if entry[0] == src_pos:
                        break
                    remove_index += 1
                candidates[id2].pop(remove_index)
                # print(f"after removing {src_pos} cand_pos is\n{cand_pos[id2]}\nand candidates is\n{candidates[id2]}")

            if opening:
                # Source is already settled, we need to open the correct ports to allow routing, keep the embedding AND leave space for following edges
                next_edge_pos, prev_edge_pos, next_off, prev_off = cls.next_and_prev_positions(id1, id2, map_circ_ord, map)
                grid.open_ports(src_pos, next_edge_pos, prev_edge_pos, next_off, prev_off)
            
            ### ASK WHY IN THE CODE BELOW THIS CONDITION GOES FIRST? DOES IT MATTERS?
            
            # Setting offset weights for candidate target positions
            tgt_pos = OctiMetroHeuristicSchematic.prepare_unsettled_node(candidates, id2, grid, cand_positions2)
            added_tgt = True

        elif src_pos == -1:
            # print("source not settled")
            # removing tgt, from src candidates
            # cand_positions1 = cand_pos[id1]
            cand_positions1 = cand_pos_src
            if tgt_pos in cand_positions1:
                cand_positions1.remove(tgt_pos)

            # Setting offset weights for candidate source positions
            src_pos = OctiMetroHeuristicSchematic.prepare_unsettled_node(candidates, id1, grid, cand_positions1)
            added_src = True

            if opening:
                print()
                #Check the neighbors of the source position
                #cls.draw_network_state(grid.graph,[tgt_pos])
                # Target is already settled, we need to open the correct ports to allow routing, keep the embedding AND leave space for following edges
                next_edge_pos, prev_edge_pos, next_off, prev_off = cls.next_and_prev_positions(id2, id1, map_circ_ord, map)
                grid.open_ports(tgt_pos, next_edge_pos, prev_edge_pos, next_off, prev_off)
                ### ADDED FOR DEBUGGING PURPOSES ###########
                #cls.draw_network_state(grid.graph,[tgt_pos])

        else:
            #### ADDED FOR DEBBUGGING PURPOSES #############
            #cls.draw_network_state(grid.graph,[src_pos,tgt_pos])

            # print("both settled")
            if opening:
                if map.graph.degree(id1) == 1:
                    grid.open_all_ports(src_pos)
                else:
                    # Source is already settled, we need to open the correct ports to allow routing, keep the embedding AND leave space for following edges
                    # print("nextingA")
                    next_edge_pos, prev_edge_pos, next_off, prev_off = cls.next_and_prev_positions(id1, id2, map_circ_ord, map)
                    # print("openingA")
                    grid.open_ports(src_pos, next_edge_pos, prev_edge_pos, next_off, prev_off)

            if opening:
                if map.graph.degree(id2) == 1:
                    grid.open_all_ports(tgt_pos)
                else:
                    # Target is already settled, we need to open the correct ports to allow routing, keep the embedding AND leave space for following edges
                    # print("nextingB")
                    next_edge_pos, prev_edge_pos, next_off, prev_off = cls.next_and_prev_positions(id2, id1, map_circ_ord, map)
                    # print("openingB")
                    grid.open_ports(tgt_pos, next_edge_pos, prev_edge_pos, next_off, prev_off)

            #### ADDED FOR DEBBUGGING PURPOSES #############
            #cls.draw_network_state(grid.graph,[src_pos,tgt_pos])

            # print("done")
            # length, path = nx.single_source_dijkstra(grid.graph, src_pos, tgt_pos)

        # Regardless of added temporary nodes, we now only need to find a shortest path between src_pos and tgt_pos
        try:
            length, path = nx.single_source_dijkstra(grid.graph, src_pos, tgt_pos)
        except nx.exception.NetworkXNoPath:
            print("NetworkX Error")
            return []

        ########## ADDED TO UNDERSTAND HOW THE PATH IS CREATED BETWEEN SOURCE AND TARGET NODES ############
        
        temp_g=[]
        for node in path:
            temp_g.append(node)
            neighbors=[n for n in grid.graph.neighbors(node)]    
            temp_g=temp_g+neighbors
        
        temp_g=list(set(temp_g))

        H = nx.subgraph(grid.graph, temp_g)
         
        pos = {}
        for node, data in H.nodes(data=True):
            # print(f"DEBUG: {node} has {data}")
            pos[node] = np.array([data["x"], data["y"]])

        labels = nx.get_edge_attributes(H,'weight')
        for label,value in labels.items():
            labels[label]=np.round(value,3)

        colors = [ 'red' if n in path else 'lightblue'  for n, data in H.nodes(data=True) ]

        #nx.draw(H,pos, with_labels=True, font_weight='bold', node_color=colors, node_size=200,connectionstyle='arc3, rad = 0.1')
        #nx.draw_networkx_edge_labels(H,pos,edge_labels=labels)

        
        ##################################################################################################


        # We now need to remove the temporary node both from the path and the network again IF we added them
        if added_src:
            grid.remove_star(src_pos)
            path = path[1:]
        if added_tgt:
            grid.remove_star(tgt_pos)
            path = path[:-1]

        ########## ADDED TO UNDERSTAND WHAT HAPPENS TO THE THE PATH AFTER REMOVING TEMPORARY NODES ############
        
        temp_g=[]
        for node in path:
            temp_g.append(node)
            neighbors=[n for n in grid.graph.neighbors(node)]    
            temp_g=temp_g+neighbors
        
        temp_g=list(set(temp_g))

        H = nx.subgraph(grid.graph, temp_g)
         
        pos = {}
        for node, data in H.nodes(data=True):
            # print(f"DEBUG: {node} has {data}")
            pos[node] = np.array([data["x"], data["y"]])

        labels = nx.get_edge_attributes(H,'weight')
        for label,value in labels.items():
            labels[label]=np.round(value,3)

        colors = [ 'red' if n in path else 'lightblue'  for n, data in H.nodes(data=True) ]

        #nx.draw(H,pos, with_labels=True, font_weight='bold', node_color=colors, node_size=200,connectionstyle='arc3, rad = 0.1')
        #nx.draw_networkx_edge_labels(H,pos,edge_labels=labels)

        
        ##################################################################################################

        return path

    @classmethod
    def remove_candidate_overlap(cls, list1, list2, pos_1, pos_2, grid):
        '''
        Removes the overlap in two provided candidate sets based on euclidean diswtance
        '''
        set1, set2 = set(list1), set(list2)
        remove1, remove2 = set(), set()
        for pos in set1:
            if pos in set2:
                posx, posy = grid.nodes[pos]["x"], grid.nodes[pos]["y"]
                if dist(pos_1[0], pos_1[1], posx, posy) < dist(pos_2[0], pos_2[1], posx, posy):
                    remove2.add(pos)
                else:
                    remove1.add(pos)
        return1, return2 = [], []
        for pos in list1:
            if pos not in remove1:
                return1.append(pos)
        for pos in list2:
            if pos not in remove2:
                return2.append(pos)
        return return1, return2

    @classmethod
    def compute_cand_poss(cls, map, options):
        # We compute the candidate positions for all stations in the new grid Since we are using the 'distances'
        # flag, we also obtain the distance of every candidate to the station as the second entry in every tuple
        cand_set_options = {
            "grid": options["grid"],
            "k": options["cand_set_factor"] if "cand_set_factor" in options else 1,
            "method": options["cand_set_method"] if "cand_set_method" in options else "circle",
            "cell_size": options["cell_size"],
            "distances": True
        }
        if "candidate_ignore" in options:
            cand_set_options["ignore"] = options["candidate_ignore"]
        if "prefer_as_candidate" in options:
            cand_set_options["prefer_attribute"] = options["prefer_as_candidate"]
            cand_set_options["prefer_factor"] = options["prefer_as_candidate_factor"]
        candidates = GridPlacementUtility.execute(map, cand_set_options)
        cand_pos = {}
        for node in candidates:
            smooth_adj = False
            for neighbor in map.graph.neighbors(node):
                edge = map.graph[node][neighbor]["edge"]
                # didnt understand quite well
                if is_smooth(edge):
                    smooth_adj = True

            candidates[node] = sorted(candidates[node], key=lambda x: x[1])

            # Kinda specific. TODO: refactor this somehow into the gridplacementutility as some kind of:
            #  TODO: if you have this contextual information as a node ("smooth"), you only get this many candidates
            if smooth_adj:
                # print(f"{node} is smooth")
                cand_pos[node] = [i[0] for i in candidates[node][:options["prefer_cand_count"] if "prefer_cand_count" in options else 1]]
                candidates[node] = candidates[node][:options["prefer_cand_count"] if "prefer_cand_count" in options else 1]
            else:
                cand_pos[node] = [i[0] for i in candidates[node]]
        return candidates

    # TODO: THIS ENTIRE FUNCTION NEEDS TO BE DEBUGGED AND VERIFIED
    @classmethod
    def perform_local_search(cls, map, grid, simple_grid, candidates, map_circ_ord, attempts=0):
        print(f"settled_edges:\n{cls.settled_edges}")

        # render_grid(grid.graph, f"{cls.OUT_DIR}/8_{cls.NAME}_grid_LS.png", dpi=600)

        counter = 1
        print(
            f"Calling local search with\nsettled_edges: {cls.settled_edges.keys()}\nsettled_edge_costs: {cls.settled_edges_costs}\nsettled_nodes: {cls.settled_nodes}")

        # print("grid is:")
        # for u, v, data in grid.graph.edges(data=True):
        #     print(data["weight"])
        # print()

        change = True
        while change:
            change = False

            # breaking = input("Press x to terminate local search\n")
            # if breaking == 'x':
            #     break

            counting = 0

            for node in map.graph.nodes():
                print(f"redoing {node}")

                if node not in cls.settled_nodes:
                    print("skipping an unsettled node (when does this happen?)")
                    continue

                if is_smooth_station(map.getStation(node)):
                    print(f"skipping smooth station {node}")
                    continue

                # nghbs = map.graph.neighbors(node)
                # Remember status quo

                grid_nghbs = simple_grid.neighbors(cls.settled_nodes[node])

                for gn in grid_nghbs:

                    # if gn not in [x[0] for x in candidates[node]]:
                    #     print("not moving out of candidates")
                    #     continue

                    if grid.graph.nodes()[gn]["shape_part"]:
                        print("Not moving onto shape")
                        continue

                    counting += 1
                    duplicate = False
                    for n in map.graph.nodes():
                        if n in cls.settled_nodes:
                            if cls.settled_nodes[n] == gn:
                                print(f"grid position is already occupied by {n}")
                                duplicate = True
                    # print(f"trying {node} at {gn}")

                    if duplicate:
                        continue

                    old_cost = 0
                    for n in map.graph.neighbors(node):
                        if (node, n) in cls.settled_edges or (n, node) in cls.settled_edges:
                            old_cost += cls.get_from_settled_edge_cost(node, n)
                                # cls.settled_edges_costs[(node, n)]
                        # elif (n, node) in cls.settled_edges:
                        #     old_cost += cls.settled_edges_costs[(n, node)]
                        else:
                            old_cost = math.inf

                    # print(f"\told cost: {old_cost}")

                    old_node_pos = cls.settled_nodes[node]
                    # old_grid = copy.deepcopy(grid)
                    # old_settled_edges = copy.deepcopy(cls.settled_edges)
                    # old_settled_edges_costs = copy.deepcopy(cls.settled_edges_costs)
                    # old_map = copy.deepcopy(map)

                    # ports_to_close = set()
                    # ports_to_open = set()
                    paths_to_close = {}
                    paths_to_open = {}
                    control_points_to_reset = {}
                    # nodes_to_unsettle = set()


                    failure = False

                    new_cost = 0
                    cls.settled_nodes[node] = gn

                    # MatplotlibGridRender.render(map, {
                    #     "grid": simple_grid,
                    #     "just_grid": True,
                    #     "grid_labels": True,
                    #     # "grid_edge_labels": "weight",
                    #     "show_nodes": True,
                    #     "font_size": 5,
                    #     # "highlight_nodes": high_nodes,
                    #     "dpi": 450,
                    #     "to_file": f"{cls.OUT_DIR}/LS/{counting}_start.png"
                    # })

                    # print(f"opening connections for {node} at {old_node_pos}")
                    grid.open_node(old_node_pos)
                    # remove previous routing of all adjacent edges, if present
                    for n in map.graph.neighbors(node):

                        if is_smooth(map.getEdge(node, n)):
                            print("skipping smooth edge")
                            continue
                        # if node == map.getEdge(node, n).source:
                        #     id1 = node
                        #     id2 = n
                        # else:
                        #     id1 = n
                        #     id2 = node

                        # if n in cls.settled_nodes:
                            # print(f"\topening connections for {n} at {cls.settled_nodes[n]}")
                        old_path = []
                        if (node, n) in cls.settled_edges or (n, node) in cls.settled_edges:
                            old_path = cls.get_from_settled_edges(node, n)
                            cls.remove_from_settled_edges(node, n)
                                # cls.settled_edges.pop((id1, id2))
                            old_path_cost = cls.get_from_settled_edge_cost(node, n)
                            cls.remove_from_settled_edge_cost(node, n)
                                # cls.settled_edges_costs.pop((id1, id2))
                        paths_to_close[(node, n)] = (old_path, old_path_cost)   # This is to reset properly
                        control_points_to_reset[(node, n)] = map.graph[node][n]["edge"].controlPoints

                        if len(old_path) > 0:
                            for i in range(len(old_path[:-1])):
                                grid.open_node(old_path[i + 1])
                                grid.open_edge(old_path[i], old_path[i + 1])
                            grid.open_edge(old_path[-2], old_path[-1])

                    # important to skip opening the ports
                    second = False
                    debug_string = "\n"
                    for n in map.graph.neighbors(node):

                        # if node == map.getEdge(node, n).source:
                        #     id1 = node
                        #     id2 = n
                        # else:
                        #     id1 = n
                        #     id2 = node

                        # value1, value2 = int(node[1:]), int(n[1:])
                        # firstID = "n" + str(min(value1, value2))
                        # secondID = "n" + str(max(value1, value2))

                        print(f"\trerouting {(node, n)}")
                        path = cls.route_single_edge(map, grid, simple_grid, node, n, map_circ_ord, candidates, second)
                        # print(f"\trouting done")
                        second = True
                        if len(path) == 0:
                            new_cost = math.inf
                            failure = True
                            # print(f"\tcould not route {(node, n)}")
                            break
                        path_cost = 0
                        for i in range(1, len(path)):
                            path_cost += grid.graph[path[i - 1]][path[i]]["weight"]
                        # print(f"\trerouted {(node, n)} for {path_cost}")

                        # print(f"\tNew Path:")
                        debug_string += f"{(node, n)}:\n"
                        for pathnode in path:
                            if not grid.is_port(pathnode):
                                # bprint(pathnode, end=", ")
                                debug_string += f"{pathnode}, "
                            else:
                                # print(f"{pathnode} ({grid.get_parent(pathnode)})", end=", ")
                                debug_string += f"{pathnode} ({grid.get_parent(pathnode)}), "
                        debug_string += "\n"

                        new_cost += path_cost

                        paths_to_open[n] = path

                        # Keep track of settled nodes and edges
                        cls.settled_nodes[node] = path[0]
                        cls.settled_nodes[n] = path[-1]

                        cls.set_settled_edge(node, n, path)
                        # cls.settled_edges[(id1, id2)] = path
                        # cls.settled_edges[(n, node)] = path[::-1]

                        cls.set_settled_edge_cost(node, n, path_cost)
                        # cls.settled_edges_costs[(id1, id2)] = path_cost
                        # cls.settled_edges_costs[(n, node)] = path_cost

                        # update station positions
                        station1, station2 = map.graph.nodes[node]["station"], map.graph.nodes[n]["station"]
                        station1.x, station1.y = grid.coord(path[0])
                        station2.x, station2.y = grid.coord(path[-1])

                        control_points = []
                        for pathnode in path[2:-2:2]:
                            control_points.append(grid.coord(pathnode))
                        map.getEdge(node, n).setControlPoints(control_points, start=node)

                        # map.graph[node][n]["edge"].controlPoints = control_points

                        # Close edges and nodes, which were used by this path
                        grid.close_path(path)
                        grid.close_node(path[0])
                        grid.close_node(path[-1])

                    # print(f"\told cost: {old_cost}")
                    # print(f"\tnewcost: {new_cost}")

                    if failure or math.isinf(new_cost) or new_cost >= old_cost:
                        # resetting to old values
                        # print("\tNo improvement")

                        # cls.settled_edges = old_settled_edges
                        # cls.settled_edges_costs = old_settled_edges_costs
                        cls.settled_nodes[node] = old_node_pos

                        # map = old_map
                        for n in map.graph.neighbors(node):
                            map.graph[node][n]["edge"].controlPoints = control_points_to_reset[(node, n)]
                        station_to_reset = map.graph.nodes[node]["station"]
                        station_to_reset.x, station_to_reset.y = grid.coord(old_node_pos)

                        grid.open_node(gn)
                        for neighbor in paths_to_open:
                            path = paths_to_open[neighbor]
                            grid.open_path(path)
                            # if path[0] != gn:
                            #     grid.open_node(path[0])
                            #     grid.open_edge(path[1], path[2])

                        for edge in paths_to_close:
                            previous_path, previous_cost = paths_to_close[edge][0], paths_to_close[edge][1]

                            cls.set_settled_edge(edge[0], edge[1], previous_path)
                            # cls.settled_edges[edge] = previous_path
                            cls.set_settled_edge_cost(edge[0], edge[1], previous_cost)
                            # cls.settled_edges_costs[edge] = previous_cost

                            # cls.settled_edges[(edge[1], edge[0])] = previous_path[::-1]
                            # cls.settled_edges_costs[(edge[1], edge[0])] = previous_cost

                            grid.close_path(previous_path)
                            grid.close_node(previous_path[0])
                            grid.close_node(previous_path[-1])

                        # grid = old_grid

                        second = False
                    else:
                        change = True
                        print(f"\tMoved {node} from {old_node_pos} to {cls.settled_nodes[node]}")
                        print(f"paths are: {debug_string}")
                        render_map(map, grid.graph, f"{cls.OUT_DIR}/8_{cls.NAME}_new_LS_{counter}.png")
                        counter += 1
                        # input()

    @classmethod
    def perform_local_search_alt(cls, map, grid, simple_grid, candidates, map_circ_ord, attempts=0):
        
        print(f"settled_edges:\n{cls.settled_edges}")

        counter = 1
        print(
            f"Calling local search with\nsettled_edges: {cls.settled_edges.keys()}\nsettled_edge_costs: {cls.settled_edges_costs}\nsettled_nodes: {cls.settled_nodes}")


        change = True
        while change:
            change = False

            counting = 0

            for node in map.graph.nodes():

                print(f"redoing {node}")

                if node not in cls.settled_nodes:
                    print("skipping an unsettled node (when does this happen?)")
                    continue

                if is_smooth_station(map.getStation(node)):
                    print(f"skipping smooth station {node}")
                    continue

                # Draw network state of current node 
                cls.draw_network_state(grid.graph,[cls.settled_nodes[node]])

                grid_nghbs = simple_grid.neighbors(cls.settled_nodes[node])

                for gn in grid_nghbs:

                    if grid.graph.nodes()[gn]["shape_part"]:
                        print("Not moving onto shape")
                        continue

                    counting += 1
                    duplicate = False

                    ### PROPOSAL ####
                    if gn in cls.settled_nodes.values():
                        n= list(cls.settled_nodes.keys())[list(cls.settled_nodes.values()).index(gn)] 
                        print("grid position is already occupied by {}".format(n))
                        duplicate=True

                    if duplicate:
                        continue

                    ### PROPOSAL (EXTRACT WEIGHTS BEFORE MODIFYING ANYTHING) ####
                    original_weights = nx.get_edge_attributes(grid.graph,'weight')


                    ### JUST FOR DEBUGGING PURPOSES ###
                    # original_weights_dict={}
                    # for key,value in original_weights.items():
                    #     original_weights_dict[','.join(str(key))]=value

                    # with open('convert.txt', 'w') as convert_file:
                    #     convert_file.write(json.dumps(original_weights_dict))

                    #### SHOULD BE IMPROVE TIME NO NEED TO CALCULATE AGAIN IF LAST TIME THERE WAS A FAILURE #####
                                        
                    old_cost = 0
                    
                    
                    for n in map.graph.neighbors(node):
                        if (node, n) in cls.settled_edges or (n, node) in cls.settled_edges:
                            old_cost += cls.get_from_settled_edge_cost(node, n)
                        else:
                            old_cost = math.inf

                    
                    old_node_pos = cls.settled_nodes[node]

                    paths_to_close = {}
                    paths_to_open = {}
                    control_points_to_reset = {}

                    failure = False

                    new_cost = 0
                    cls.settled_nodes[node] = gn

                    grid.open_node(old_node_pos)
                    # remove previous routing of all adjacent edges, if present

                    cls.draw_network_state(grid.graph,[old_node_pos])

                    for n in map.graph.neighbors(node):

                        if is_smooth(map.getEdge(node, n)):
                            print("skipping smooth edge")
                            continue

                        old_path = []
                        if (node, n) in cls.settled_edges or (n, node) in cls.settled_edges:
                            old_path = cls.get_from_settled_edges(node, n)
                            cls.remove_from_settled_edges(node, n)
                            old_path_cost = cls.get_from_settled_edge_cost(node, n)
                            cls.remove_from_settled_edge_cost(node, n)

                        paths_to_close[(node, n)] = (old_path, old_path_cost)   # This is to reset properly
                        control_points_to_reset[(node, n)] = map.graph[node][n]["edge"].controlPoints

                        cls.draw_network_state(grid.graph,old_path)


                        if len(old_path) > 0:
                            for i in range(len(old_path[:-1])):
                                grid.open_node(old_path[i + 1])
                                grid.open_edge(old_path[i], old_path[i + 1])
                                cls.draw_network_state(grid.graph,[old_path[i+1]])
                            grid.open_edge(old_path[-2], old_path[-1])

                        cls.draw_network_state(grid.graph,old_path)

                    # important to skip opening the ports
                    second = False
                    debug_string = "\n"
                    for n in map.graph.neighbors(node):


                        print(f"\trerouting {(node, n)}")
                        path = cls.route_single_edge(map, grid, simple_grid, node, n, map_circ_ord, candidates, second)
                        # print(f"\trouting done")
                        
                        second = True
                        if len(path) == 0:
                            new_cost = math.inf
                            failure = True
                            # print(f"\tcould not route {(node, n)}")
                            break
                        path_cost = 0
                        for i in range(1, len(path)):
                            path_cost += grid.graph[path[i - 1]][path[i]]["weight"]
                        # print(f"\trerouted {(node, n)} for {path_cost}")

                        # print(f"\tNew Path:")
                        debug_string += f"{(node, n)}:\n"
                        for pathnode in path:
                            if not grid.is_port(pathnode):
                                # bprint(pathnode, end=", ")
                                debug_string += f"{pathnode}, "
                            else:
                                # print(f"{pathnode} ({grid.get_parent(pathnode)})", end=", ")
                                debug_string += f"{pathnode} ({grid.get_parent(pathnode)}), "
                        debug_string += "\n"

                        new_cost += path_cost

                        paths_to_open[n] = path

                        ###CHECK HOW ARE THE NODES OF THE FOUND PATH BEFORE ANY OPERATION TO KNOW IF AT THE END THEY ARE THE SAME
                        cls.draw_network_state(grid.graph,path)


                        # Keep track of settled nodes and edges
                        cls.settled_nodes[node] = path[0]
                        cls.settled_nodes[n] = path[-1] 

                        cls.set_settled_edge(node, n, path)
                        # cls.settled_edges[(id1, id2)] = path
                        # cls.settled_edges[(n, node)] = path[::-1]

                        cls.set_settled_edge_cost(node, n, path_cost)
                        # cls.settled_edges_costs[(id1, id2)] = path_cost
                        # cls.settled_edges_costs[(n, node)] = path_cost

                        # update station positions
                        station1, station2 = map.graph.nodes[node]["station"], map.graph.nodes[n]["station"]
                        station1.x, station1.y = grid.coord(path[0])
                        station2.x, station2.y = grid.coord(path[-1])

                        control_points = []
                        for pathnode in path[2:-2:2]:
                            control_points.append(grid.coord(pathnode))
                        map.getEdge(node, n).setControlPoints(control_points, start=node)

                        # Close edges and nodes, which were used by this path
                        grid.close_path(path)
                        grid.close_node(path[0])
                        grid.close_node(path[-1])

                        ###CHECK HOW ARE THE NODES OF THE FOUND PATH AFTER OPERATION TO KNOW IF AT THE END THEY ARE THE SAME
                        cls.draw_network_state(grid.graph,path)

                        ###ADDED TO AVOID UNNECESARY LOOPS
                        if new_cost == math.inf:
                           break
                    
                    if failure or math.isinf(new_cost) or new_cost >= old_cost:

                        nx.set_edge_attributes(grid.graph,original_weights,name='weight')

                        # modified_weights = nx.get_edge_attributes(grid.graph,'weight')

                        # JUST FOR DEBUGGING PURPOSES
                        # modified_weights_dict={}
                        # for key,value in modified_weights.items():
                        #         modified_weights_dict[','.join(str(key))]=value

                        # with open('convert2.txt', 'w') as convert_file:
                        #     convert_file.write(json.dumps(modified_weights_dict))


                        cls.settled_nodes[node] = old_node_pos

                        
                        for n in map.graph.neighbors(node):
                            map.graph[node][n]["edge"].controlPoints = control_points_to_reset[(node, n)]
                        station_to_reset = map.graph.nodes[node]["station"]
                        station_to_reset.x, station_to_reset.y = grid.coord(old_node_pos)

                        #grid.open_node(gn)
                        for neighbor in paths_to_open:
                            #path = paths_to_open[neighbor]
                            #grid.open_path(path)

                            ## CHECK MODIFIED PATH IS SAME AS BEFORE
                            cls.draw_network_state(grid.graph,path)


                            # if path[0] != gn:
                            #     grid.open_node(path[0])
                            #     grid.open_edge(path[1], path[2])

                        


                        for edge in paths_to_close:
                            previous_path, previous_cost = paths_to_close[edge][0], paths_to_close[edge][1]

                            cls.set_settled_edge(edge[0], edge[1], previous_path)
                            cls.set_settled_edge_cost(edge[0], edge[1], previous_cost)
                            
                            # cls.settled_edges_costs[edge] = previous_cost

                            # cls.settled_edges[(edge[1], edge[0])] = previous_path[::-1]
                            # cls.settled_edges_costs[(edge[1], edge[0])] = previous_cost

                        #     grid.close_path(previous_path)
                        #     grid.close_node(previous_path[0])
                        #     grid.close_node(previous_path[-1])

                        ## CHECK STATE OF ORIGINAL NODE IS SAME AS BEFORE
                        cls.draw_network_state(grid.graph,[old_node_pos])

                        cls.draw_network_state(grid.graph,old_path)


                    #     # grid = old_grid

                        second = False
                    else:
                        change = True
                        print(f"\tMoved {node} from {old_node_pos} to {cls.settled_nodes[node]}")
                        print(f"paths are: {debug_string}")
                        render_map(map, grid.graph, f"{cls.OUT_DIR}/8_{cls.NAME}_new_LS_{counter}.png")
                        counter += 1
                        # input()


    # TODO: THIS ENTIRE FUNCTION NEEDS TO BE DEBUGGED AND VERIFIED
    @classmethod
    def perform_local_search_dc(cls, map, grid, simple_grid, candidates, map_circ_ord, attempts=0):
        counter = 1
        print(
            f"Calling local search with\nsettled_edges: {cls.settled_edges.keys()}\nsettled_edge_costs: {cls.settled_edges_costs}\nsettled_nodes: {cls.settled_nodes}")

        # print("grid is:")
        # for u, v, data in grid.graph.edges(data=True):
        #     print(data["weight"])
        # print()

        change = True
        while change:
            change = False

            breaking = input("Press x to terminate local search\n")
            if breaking == 'x':
                break

            counting = 0

            for node in map.graph.nodes():
                print(f"redoing {node}")

                if node not in cls.settled_nodes:
                    print("skipping an unsettled node (when does this happen?)")
                    continue

                if is_smooth_station(map.getStation(node)):
                    print(f"skipping smooth station {node}")
                    continue

                # nghbs = map.graph.neighbors(node)
                # Remember status quo

                grid_nghbs = simple_grid.neighbors(cls.settled_nodes[node])

                for gn in grid_nghbs:
                    counting += 1
                    duplicate = False
                    for n in map.graph.nodes():
                        if n in cls.settled_nodes:
                            if cls.settled_nodes[n] == gn:
                                print(f"grid position is already occupied by {n}")
                                duplicate = True
                    # print(f"trying {node} at {gn}")

                    if duplicate:
                        continue

                    old_cost = 0
                    for n in map.graph.neighbors(node):
                        if (node, n) in cls.settled_edges or (n, node) in cls.settled_edges:
                            old_cost += cls.get_from_settled_edge_cost(node, n)
                                # cls.settled_edges_costs[(node, n)]
                        # elif (n, node) in cls.settled_edges:
                        #     old_cost += cls.settled_edges_costs[(n, node)]
                        else:
                            old_cost = math.inf

                    # print(f"\told cost: {old_cost}")

                    old_node_pos = cls.settled_nodes[node]
                    old_grid = copy.deepcopy(grid)
                    old_settled_edges = copy.deepcopy(cls.settled_edges)
                    old_settled_edges_costs = copy.deepcopy(cls.settled_edges_costs)
                    old_map = copy.deepcopy(map)

                    failure = False

                    new_cost = 0
                    cls.settled_nodes[node] = gn

                    # print(f"opening connections for {node} at {old_node_pos}")
                    grid.open_node(old_node_pos)
                    # remove previous routing of all adjacent edges, if present
                    for n in map.graph.neighbors(node):
                        # if n in cls.settled_nodes:
                        # print(f"\topening connections for {n} at {cls.settled_nodes[n]}")
                        old_path = []
                        if (node, n) in cls.settled_edges or (n, node) in cls.settled_edges:
                            old_path = cls.get_from_settled_edges(node, n)
                            cls.remove_from_settled_edges(node, n)
                                # cls.settled_edges.pop((id1, id2))
                            old_path_cost = cls.get_from_settled_edge_cost(node, n)
                            cls.remove_from_settled_edge_cost(node, n)

                        # if (n, node) in cls.settled_edges:
                        #     old_path = cls.settled_edges.pop((n, node))
                        #     old_path_cost = cls.settled_edges_costs.pop((n, node))

                        if len(old_path) > 0:
                            for i in range(len(old_path[1:-1])):
                                grid.open_node(old_path[i + 1])
                                grid.open_edge(old_path[i], old_path[i + 1])
                            grid.open_edge(old_path[-2], old_path[-1])

                    # important to skip opening the ports
                    second = False
                    debug_string = "\n"
                    for n in map.graph.neighbors(node):
                        print(f"\trerouting {(node, n)}")
                        path = cls.route_single_edge(map, grid, simple_grid, node, n, map_circ_ord, candidates,
                                                     second)
                        # print(f"\trouting done")
                        second = True
                        if len(path) == 0:
                            new_cost = math.inf
                            failure = True
                            # print(f"\tcould not route {(node, n)}")
                            break
                        path_cost = 0
                        for i in range(1, len(path)):
                            path_cost += grid.graph[path[i - 1]][path[i]]["weight"]
                        # print(f"\trerouted {(node, n)} for {path_cost}")

                        # print(f"\tNew Path:")
                        debug_string += f"{(node, n)}:\n"
                        for pathnode in path:
                            if not grid.is_port(pathnode):
                                # bprint(pathnode, end=", ")
                                debug_string += f"{pathnode}, "
                            else:
                                # print(f"{pathnode} ({grid.get_parent(pathnode)})", end=", ")
                                debug_string += f"{pathnode} ({grid.get_parent(pathnode)}), "
                        debug_string += "\n"

                        new_cost += path_cost

                        # Keep track of settled nodes and edges
                        cls.settled_nodes[node] = path[0]
                        cls.settled_nodes[n] = path[-1]

                        cls.set_settled_edge(node, n, path)
                        # cls.settled_edges[(node, n)] = path
                        # cls.settled_edges[(n, node)] = path[::-1]
                        cls.set_settled_edge_cost(node, n, path_cost)
                        # cls.settled_edges_costs[(node, n)] = path_cost
                        # cls.settled_edges_costs[(n, node)] = path_cost

                        # update station positions
                        station1, station2 = map.graph.nodes[node]["station"], map.graph.nodes[n]["station"]
                        station1.x, station1.y = grid.coord(path[0])
                        station2.x, station2.y = grid.coord(path[-1])

                        control_points = []
                        for pathnode in path[2:-2:2]:
                            control_points.append(grid.coord(pathnode))


                        map.getEdge(node, n).setControlPoints(control_points, start=node)

                        map.graph[node][n]["edge"].controlPoints = control_points

                        # Close edges and nodes, which were used by this path
                        grid.close_path(path)
                        grid.close_node(path[0])
                        grid.close_node(path[-1])

                    # print(f"\told cost: {old_cost}")
                    # print(f"\tnewcost: {new_cost}")

                    if failure or math.isinf(new_cost) or new_cost >= old_cost:
                        # resetting to old values
                        # print("\tNo improvement")

                        cls.settled_edges = old_settled_edges
                        cls.settled_edges_costs = old_settled_edges_costs

                        map = old_map
                        grid = old_grid

                        second = False
                    else:
                        change = True
                        print(f"\tMoved {node} from {old_node_pos} to {cls.settled_nodes[node]}")
                        print(f"paths are: {debug_string}")
                        render_map(map, grid.graph, f"{cls.OUT_DIR}/8_{cls.NAME}_new_LS_{counter}.png")
                        counter += 1
                        # input()


    # TODO: the three methods below are technically grid methods? IF so move to extended grid class
    @classmethod
    def prepare_unsettled_node(cls, candidates, id1, grid, cand_positions1, original_distance_factor=1):
        # Setting offset weights for candidate source positions
        weights = {}
        for candidate in candidates[id1]:
            weights[candidate[0]] = candidate[1] * original_distance_factor
        # Adding connections to a one off source node accordingly
        return grid.add_star(cand_positions1, weights, cls.orig_positions[id1])

    @classmethod
    def next_and_prev_positions(cls, id1, id2, map_circ_ord, map):
        '''
        Compute the next and the previous edge, which are already settled. This edge has to exist. Can be the same.
        '''
        # print(f"getting next and previous for {id1} and {id2}")
        # print("next?")
        next_set_target, next_off = cls.next_settled(id1, id2, map_circ_ord)
        # print(f"next settled target is {next_set_target} with path {cls.get_from_settled_edges(id1, next_set_target)}")
        # print("previous?")
        prev_set_target, prev_off = cls.next_settled(id1, id2, map_circ_ord, backwards=True)
        # print(f"previous settled target is {prev_set_target} with path {cls.get_from_settled_edges(id1, prev_set_target)}")

        # print(f"id1: {id1} next_set_target: {next_set_target} prev_set_target: {prev_set_target}")

        value1, value2, value3 = int(id1[1:]), int(next_set_target[1:]), int(prev_set_target[1:])

        # print(f"v1: {value1}, v2: {value2}, v3: {value3}")

        # next_edge_pos = cls.get_from_settled_edges(id1, next_set_target)[1]
        # prev_edge_pos = cls.get_from_settled_edges(id1, prev_set_target)[1]
        next_path = cls.get_from_settled_edges(id1, next_set_target)
        id1_pos = cls.settled_nodes[id1]
        if id1_pos == next_path[0]:
            next_edge_pos = next_path[1]
        elif id1_pos == next_path[-1]:
            next_edge_pos = next_path[-2]
        else:
            print("This does not work")
            quit()

        # if value1 < value2:
        #     next_edge_pos = cls.get_from_settled_edges(id1, next_set_target)[1]
        # else:
        #     next_edge_pos = cls.get_from_settled_edges(id1, next_set_target)[-2]

        # cls.settled_edges[(id1, next_set_target)][1]

        prev_path = cls.get_from_settled_edges(id1, prev_set_target)
        if id1_pos == prev_path[0]:
            prev_edge_pos = prev_path[1]
        elif id1_pos == prev_path[-1]:
            prev_edge_pos = prev_path[-2]
        else:
            print("This does not work")
            quit()

        # if value1 < value3:
        #     prev_edge_pos = cls.get_from_settled_edges(id1, prev_set_target)[1]
        # else:
        #     prev_edge_pos = cls.get_from_settled_edges(id1, prev_set_target)[-2]

        # cls.settled_edges[(id1, prev_set_target)][1]

        # if id1 == "n18":
        #     print(f"{id1} - {next_set_target} is currently routed like this: {settled_edges[(id1, next_set_target)]}")
        # next_edge_pos = map.getEdge(id1, next_set_target).getEdgeGeometry(id1)[1]
        # prev_edge_pos = map.getEdge(id1, prev_set_target).getEdgeGeometry(id1)[1]
        return next_edge_pos, prev_edge_pos, next_off, prev_off

    @classmethod
    def next_settled(cls, u, v, c_ord, backwards=False):
        '''
        Starting at one edge looks for the next edge which has been settled. Edge has to exist. Can be called backwards.
        '''
        u_neighb = c_ord[u]
        v_ind = u_neighb.index(v)
        modul = len(u_neighb)

        if not backwards:
            i = (v_ind + 1) % modul
        else:
            i = (v_ind - 1) % modul
        next_off = 0
        while True:
            if (u, u_neighb[i]) in cls.settled_edges or (u_neighb[i], u) in cls.settled_edges:
                next_e = u_neighb[i]
                break
            next_off += 1
            if not backwards:
                i += 1
            else:
                i -= 1
            i %= modul
        return next_e, next_off

    @classmethod
    def set_settled_edge(cls, id1, id2, value):
        value1, value2 = int(id1[1:]), int(id2[1:])
        firstID = "n" + str(min(value1, value2))
        secondID = "n" + str(max(value1, value2))
        # print(f"settling {(firstID, secondID)}")
        cls.settled_edges[(firstID, secondID)] = value

    @classmethod
    def get_from_settled_edges(cls, id1, id2):
        value1, value2 = int(id1[1:]), int(id2[1:])
        firstID = "n" + str(min(value1, value2))
        secondID = "n" + str(max(value1, value2))
        return cls.settled_edges[(firstID, secondID)]

    @classmethod
    def remove_from_settled_edges(cls, id1, id2):
        value1, value2 = int(id1[1:]), int(id2[1:])
        firstID = "n" + str(min(value1, value2))
        secondID = "n" + str(max(value1, value2))
        cls.settled_edges.pop((firstID, secondID))

    @classmethod
    def set_settled_edge_cost(cls, id1, id2, value):
        value1, value2 = int(id1[1:]), int(id2[1:])
        firstID = "n" + str(min(value1, value2))
        secondID = "n" + str(max(value1, value2))
        cls.settled_edges_costs[(firstID, secondID)] = value

    @classmethod
    def get_from_settled_edge_cost(cls, id1, id2):
        value1, value2 = int(id1[1:]), int(id2[1:])
        firstID = "n" + str(min(value1, value2))
        secondID = "n" + str(max(value1, value2))
        return cls.settled_edges_costs[(firstID, secondID)]

    @classmethod
    def remove_from_settled_edge_cost(cls, id1, id2):
        value1, value2 = int(id1[1:]), int(id2[1:])
        firstID = "n" + str(min(value1, value2))
        secondID = "n" + str(max(value1, value2))
        cls.settled_edges_costs.pop((firstID, secondID))

    @classmethod
    def draw_network_state(cls,graph,nodes):
        temp_g=[]
        for node in nodes:
            temp_g.append(node)
            neighbors=[n for n in graph.neighbors(node)]    
            temp_g=temp_g+neighbors
        
        temp_g=list(set(temp_g))

        H = nx.subgraph(graph, temp_g)
         
        pos = {}
        for node, data in H.nodes(data=True):
            # print(f"DEBUG: {node} has {data}")
            pos[node] = np.array([data["x"], data["y"]])

        labels = nx.get_edge_attributes(H,'weight')
        for label,value in labels.items():
            labels[label]=np.round(value,3)

        colors = [ 'red' if n in nodes else 'lightblue'  for n, data in H.nodes(data=True) ]

        nx.draw(H,pos, with_labels=True, font_weight='bold', node_color=colors, node_size=200,connectionstyle='arc3, rad = 0.1')
        nx.draw_networkx_edge_labels(H,pos,edge_labels=labels)