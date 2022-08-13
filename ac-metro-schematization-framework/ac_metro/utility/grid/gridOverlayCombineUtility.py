from ac_metro.utility.abstractUtility import AbstractUtility
from ac_metro.utility.custom.octi.helper import dist
import networkx as nx
import ground.base as g_base
import bentley_ottmann as bent_ott

from ac_metro.utility.custom.octi.poly_point_isect import isect_segments_include_segments


class GridOverlayer(AbstractUtility):

    @staticmethod
    def execute(map, options=None):

        method = options["method"] if "method" in options else "distance"
        if "grid" in options:
            grid = options["grid"]
        else:
            print(f"No grid was provided to {GridOverlayer.__name__}. Aborting")
            return
        if "shape" in options:
            shape = options["shape"]
        else:
            print(f"No shape was provided to {GridOverlayer.__name__}. Aborting")
            return

        if method == "distance":
            upper_threshhold = options["upper_threshhold"] if "upper_threshhold" in options else 1.5
            lower_threshhold = options["lower_threshhold"] if "lower_threshhold" in options else 0.5
            if "proximity_removal" in options:
                prox_removal = options["proximity_removal"]
            return distance_method(grid, shape, upper_threshhold, lower_threshhold,
                                   options["shape_connection_cost"] if "shape_connection_cost" in options else 1,
                                   options["proximity_removal"] if "proximity_removal" in options else -1)

        elif method == "crossing":
            if "add_crossing_of_shape" in options and options["add_crossing_of_shape"]:
                ADD_CROSSING_OF_SHAPE = True
            else:
                ADD_CROSSING_OF_SHAPE = False
            if "remove_original_shape" in options and options["remove_original_shape"]:
                REMOVE_ORIGINAL_SHAPE = True
            else:
                REMOVE_ORIGINAL_SHAPE = False

            if "g_left" in options:
                g_left = options["g_left"]
            else:
                print(f"No g_left was provided to {GridOverlayer.__name__}. Aborting")
                return
            if "g_bot" in options:
                g_bot = options["g_bot"]
            else:
                print(f"No g_bot was provided to {GridOverlayer.__name__}. Aborting")
                return
            if "g_width" in options:
                g_width = options["g_width"]
            else:
                print(f"No g_width was provided to {GridOverlayer.__name__}. Aborting")
                return
            if "g_height" in options:
                g_height = options["g_height"]
            else:
                print(f"No g_height was provided to {GridOverlayer.__name__}. Aborting")
                return
            if "cell_size" in options:
                cell_size = options["cell_size"]
            else:
                print(f"No cell_size was provided to {GridOverlayer.__name__}. Aborting")
                return
            return crossing_method(grid, g_left, g_bot, g_width, g_height, cell_size, shape,
                                   options["shape_cost"] if "shape_cost" in options else 1,
                                   options["shape_connection_cost"] if "shape_connection_cost" in options else 1,
                                   ADD_CROSSING_OF_SHAPE, REMOVE_ORIGINAL_SHAPE)

        else:
            print(f"Chosen method {method} for {str(GridOverlayer)} is not supported")
            return


# method where crossed edges are deleted and reconnection is done based on the distance of nodes
def distance_method(grid, shape, upper_thresshold, lower_threshhold, connection_cost, proximity_removal):
    print(f"Connecting Shape with cost {connection_cost}")

    _, _, crossed_grid_lines = compute_shape_grid_crossings_3(shape, grid)

    # If the shape crosses a grid edge, we remove the grid edge...
    for u, v in crossed_grid_lines:
        grid.remove_edge(u, v)

    nodes_to_remove = []
    if proximity_removal > 0:
        for g_node, g_data in grid.nodes(data=True):
            g_x, g_y = g_data["x"], g_data["y"]
            for s_node, s_data in shape.nodes(data=True):
                s_x, s_y = s_data["x"], s_data["y"]
                distance = dist(s_x, s_y, g_x, g_y)
                if distance < proximity_removal:
                    nodes_to_remove.append(g_node)
                    break
        for node in nodes_to_remove:
            grid.remove_node(node)


    # ...combining the two graphs...
    for node in shape.nodes():
        shape.nodes()[node]["shape_part"] = True
    for u, v in shape.edges():
        shape[u][v]["shape_part"] = True
    for node in grid.nodes():
        grid.nodes()[node]["shape_part"] = False
    for u, v in grid.edges():
        grid[u][v]["shape_part"] = False

    custom_grid = nx.compose(grid, shape)

    # ...and then reintroduce connections to the grid based on distance
    for s_node, s_data in shape.nodes(data=True):
        s_x, s_y = s_data["x"], s_data["y"]
        for g_node, g_data in grid.nodes(data=True):
            g_x, g_y = g_data["x"], g_data["y"]
            distance = dist(s_x, s_y, g_x, g_y)
            if lower_threshhold < distance < upper_thresshold:
                custom_grid.add_edge(s_node, g_node, crossing=[], shape_part=False, weight=connection_cost)
                # TODO: should shape_part be True here?

    # We mark nodes not on the shape as such
    for node, data in custom_grid.nodes(data=True):
        if "shape_part" not in data:
            data["shape_part"] = False

    # Sanity check to see if all edges have an associated weight
    for u, v, data in custom_grid.edges(data=True):
        if "weight" not in data:
            print(f"{(u, v)} has no edge data for weight")

    return custom_grid


# New crossing method
# TODO: CHECK HOW THIS WORKS
def compute_shape_grid_crossings_2(shape, grid):
    coords_to_type, coords_to_edge = {}, {}
    segments = []
    for s_u, s_v in shape.edges():
        s_u_x, s_u_y = shape.nodes()[s_u]["x"], shape.nodes()[s_u]["y"]
        s_v_x, s_v_y = shape.nodes()[s_v]["x"], shape.nodes()[s_v]["y"]

        coords_to_type[(s_u_x, s_u_y), (s_v_x, s_v_y)] = "shape"
        coords_to_edge[(s_u_x, s_u_y), (s_v_x, s_v_y)] = (s_u, s_v)
        coords_to_type[(s_v_x, s_v_y), (s_u_x, s_u_y)] = "shape"
        coords_to_edge[(s_v_x, s_v_y), (s_u_x, s_u_y)] = (s_u, s_v)
        segments.append(((s_u_x, s_u_y), (s_v_x, s_v_y)))

    for g_u, g_v in grid.edges():
        g_u_x, g_u_y = grid.nodes()[g_u]["x"], grid.nodes()[g_u]["y"]
        g_v_x, g_v_y = grid.nodes()[g_v]["x"], grid.nodes()[g_v]["y"]

        coords_to_type[(g_u_x, g_u_y), (g_v_x, g_v_y)] = "grid"
        coords_to_edge[(g_u_x, g_u_y), (g_v_x, g_v_y)] = (g_u, g_v)
        coords_to_type[(g_v_x, g_v_y), (g_u_x, g_u_y)] = "grid"
        coords_to_edge[(g_v_x, g_v_y), (g_u_x, g_u_y)] = (g_u, g_v)
        segments.append(((g_u_x, g_u_y), (g_v_x, g_v_y)))

    crossings = isect_segments_include_segments(segments)

    ID = grid.number_of_nodes() + shape.number_of_nodes()
    id_to_crossing = {}
    crossings_on_edges = {}
    crossed_grid_lines = set()
    for c in crossings:
        if(len(c[1]) > 2):
            print(f"multiple crossings in one point: {c}")
            for seg in c[1]:
                print(coords_to_edge[seg])
        seg1, seg2 = c[1]
        if coords_to_type[seg1] == coords_to_type[seg2]:
            continue

        id_to_crossing[ID] = c[0]
        u1, u2 = coords_to_edge[seg1]
        v1, v2 = coords_to_edge[seg2]
        if (u1, u2) not in crossings_on_edges:
            crossings_on_edges[(u1, u2)] = []
        if (v1, v2) not in crossings_on_edges:
            crossings_on_edges[(v1, v2)] = []
        crossings_on_edges[(u1, u2)].append(ID)
        crossings_on_edges[(v1, v2)].append(ID)

        if coords_to_type[seg1] == "grid":
            crossed_grid_lines.add(coords_to_edge[seg1])
        else:
            crossed_grid_lines.add(coords_to_edge[seg2])
        ID += 1

    return crossings_on_edges, id_to_crossing, crossed_grid_lines

# New New crossing method
# TODO: CHECK HOW THIS WORKS
def compute_shape_grid_crossings_3(shape, grid):
    coords_to_type, coords_to_edge = {}, {}
    segments = []
    for s_u, s_v in shape.edges():
        s_u_x, s_u_y = shape.nodes()[s_u]["x"], shape.nodes()[s_u]["y"]
        s_v_x, s_v_y = shape.nodes()[s_v]["x"], shape.nodes()[s_v]["y"]

        coords_to_type[(s_u_x, s_u_y), (s_v_x, s_v_y)] = "shape"
        coords_to_edge[(s_u_x, s_u_y), (s_v_x, s_v_y)] = (s_u, s_v)
        coords_to_type[(s_v_x, s_v_y), (s_u_x, s_u_y)] = "shape"
        coords_to_edge[(s_v_x, s_v_y), (s_u_x, s_u_y)] = (s_u, s_v)
        segments.append(((s_u_x, s_u_y), (s_v_x, s_v_y)))

    for g_u, g_v in grid.edges():
        g_u_x, g_u_y = grid.nodes()[g_u]["x"], grid.nodes()[g_u]["y"]
        g_v_x, g_v_y = grid.nodes()[g_v]["x"], grid.nodes()[g_v]["y"]

        coords_to_type[(g_u_x, g_u_y), (g_v_x, g_v_y)] = "grid"
        coords_to_edge[(g_u_x, g_u_y), (g_v_x, g_v_y)] = (g_u, g_v)
        coords_to_type[(g_v_x, g_v_y), (g_u_x, g_u_y)] = "grid"
        coords_to_edge[(g_v_x, g_v_y), (g_u_x, g_u_y)] = (g_u, g_v)
        segments.append(((g_u_x, g_u_y), (g_v_x, g_v_y)))

    crossings = isect_segments_include_segments(segments)

    ID = grid.number_of_nodes() + shape.number_of_nodes()
    id_to_crossing = {}
    crossings_on_edges = {}
    crossed_grid_lines = set()
    for c in crossings:

        id_to_crossing[ID] = c[0]
        crossed_segments = c[1]
        go_on = False
        for crossed_segment in crossed_segments:
            if coords_to_type[crossed_segment] != "grid":
                go_on = True
        if not go_on:
            continue
        for crossed_segment in crossed_segments:
            seg_tuple = (crossed_segment[0], crossed_segment[1])
            if seg_tuple not in crossings_on_edges:
                crossings_on_edges[seg_tuple] = []
            crossings_on_edges[seg_tuple].append(ID)
            if coords_to_type[crossed_segment] == "grid":
                crossed_grid_lines.add(coords_to_edge[crossed_segment])
        ID += 1

    return crossings_on_edges, id_to_crossing, crossed_grid_lines


# Possibly deprecated old method for crossings of shape and grid
# TODO: check if this works and should be used
def compute_shape_grid_crossings(shape, grid):
    context = g_base.get_context()
    Point, Segment = context.point_cls, context.segment_cls
    crossings_on_edges = {}
    id_to_crossing = {}
    ID = grid.number_of_nodes() + shape.number_of_nodes()
    for s_u, s_v in shape.edges():
        s_u_x, s_u_y = shape.nodes()[s_u]["x"], shape.nodes()[s_u]["y"]
        s_v_x, s_v_y = shape.nodes()[s_v]["x"], shape.nodes()[s_v]["y"]
        s_segment = Segment(Point(s_u_x, s_u_y), Point(s_v_x, s_v_y))

        for g_u, g_v in grid.edges():
            g_u_x, g_u_y = grid.nodes()[g_u]["x"], grid.nodes()[g_u]["y"]
            g_v_x, g_v_y = grid.nodes()[g_v]["x"], grid.nodes()[g_v]["y"]
            g_segment = Segment(Point(g_u_x, g_u_y), Point(g_v_x, g_v_y))

            ints = bent_ott.segments_intersections([s_segment, g_segment])
            if (0, 1) in ints:
                x_val = ints[(0, 1)][0].x
                if not isinstance(ints[(0, 1)][0].x, float):
                    x_val = ints[(0, 1)][0].x.numerator / ints[(0, 1)][0].x.denominator
                y_val = ints[(0, 1)][0].y
                if not isinstance(ints[(0, 1)][0].y, float):
                    y_val = ints[(0, 1)][0].y.numerator / ints[(0, 1)][0].y.denominator
                id_to_crossing[ID] = (x_val, y_val)
                if (s_u, s_v) in crossings_on_edges:
                    crossings_on_edges[(s_u, s_v)].append(ID)
                else:
                    crossings_on_edges[(s_u, s_v)] = [ID]
                if (g_u, g_v) in crossings_on_edges:
                    crossings_on_edges[(g_u, g_v)].append(ID)
                else:
                    crossings_on_edges[(g_u, g_v)] = [ID]
                ID += 1
    return crossings_on_edges, id_to_crossing


def remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, low, high, ID):
    if low > high:
        temp = high
        high = low
        low = temp
        print("This was the wrong way around")
    fset = (low, high)
    left_x, left_y = grid.nodes()[low]["x"], grid.nodes()[low]["y"]
    if fset not in edges_to_add_dict:
        edges_to_add_dict[fset] = []
    edges_to_add_dict[fset].append((ID, dist(left_x, left_y, nodes_to_add[ID][0], nodes_to_add[ID][1])))


def crossing_method(grid, g_left, g_bot, g_width, g_height, cell_size, shape, cost_shape_edge,
                    cost_shape_connection_edge, ADD_CROSSING_OF_SHAPE, REMOVE_ORIGINAL_SHAPE):
    print(f"Computing crossings between shape and grid")
    # crossings, IDs = compute_shape_grid_crossings(shape, grid)
    crossings, IDs, _ = compute_shape_grid_crossings_2(shape, grid)

    print(f"Adding shape nodes and edges to grid")
    for node, data in shape.nodes(data=True):
        grid.add_node(node, x=data["x"], y=data["y"], shape_part=True)
    for u, v, data in shape.edges(data=True):
        grid.add_edge(u, v, crossing=[], weight=cost_shape_edge)

    # We mark nodes not on the shape as such
    for node, data in grid.nodes(data=True):
        if "shape_part" not in data or not data["shape_part"]:
            grid.nodes()[node]["shape_part"] = False

    if ADD_CROSSING_OF_SHAPE:
        print(f"Reforming into subdivided edges")
        intersections_to_nodes(grid, crossings, IDs, shape, cost_shape_edge)

    print(f"Removing original shape nodes")
    # We remove the original points of the shape from the grid, as we do not want any degree 2 nodes in there

    if REMOVE_ORIGINAL_SHAPE:
        grid = remove_original_shape_nodes(shape, grid)

    # DEBUG Here we just print the computed grid once
    highlight_node = []
    for node in grid.nodes():
        if grid.nodes()[node]["shape_part"]:
            highlight_node.append(node)

    print(f"Octolinearizing new shape nodes")
    # We want to ensure that the newly added nodes of the shape are connected in an octilinear fashion to the
    octolinearize_shape_nodes(grid, g_left, g_bot, g_width, g_height, cell_size, cost_shape_connection_edge)

    # We mark nodes not on the shape as such
    for node, data in grid.nodes(data=True):
        if "shape_part" not in data or not data["shape_part"]:
            grid.nodes()[node]["shape_part"] = False

    # DEBUG Here we just print the computed grid once
    highlight_node = []
    for node in grid.nodes():
        if grid.nodes()[node]["shape_part"]:
            highlight_node.append(node)

    return grid


def intersections_to_nodes(grid, crossings, IDs, shape, cost_shape_edge):
    # We have to track these to not change the grid structure while iterating over it. Not needed in the second loop
    edges_to_remove = []
    edges_to_add = []
    nodes_to_add = {}

    # For every single edge of the grid...
    for u, v in grid.edges():
        # ...if it has an intersection...
        if (u, v) not in crossings:
            continue

        # ...we sort the intersections based on their distance to one endpoint...
        edge_crossings = crossings[(u, v)]
        edge_crossings.sort(key=lambda x: dist(IDs[x][0], IDs[x][1], grid.nodes()[u]["x"], grid.nodes()[u]["y"]))

        # ...add all intersections as vertices into the grid...
        for c in edge_crossings:
            if c not in grid:
                nodes_to_add[c] = (IDs[c][0], IDs[c][1])

        # ...add the edges u - int[0] - int[1] -...- int[len(int)-1] - v
        for i in range(len(edge_crossings) - 1):
            edges_to_add.append((edge_crossings[i], edge_crossings[i + 1]))
        edges_to_add.append((u, edge_crossings[0]))
        edges_to_add.append((edge_crossings[len(edge_crossings) - 1], v))

        # ...and finally remove the original edge from the grid
        edges_to_remove.append((u, v))

    # Actual adding of vertices and edges happens here
    for n in nodes_to_add:
        # We mark newly added nodes as part of the shape
        grid.add_node(n, x=nodes_to_add[n][0], y=nodes_to_add[n][1], shape_part=True)
    for e in edges_to_remove:
        if e[0] not in grid:
            print(f"{e[0]} is not in the grid")
        if e[1] not in grid:
            print(f"{e[1]} is not in the grid")
        grid.remove_edge(e[0], e[1])
    for e in edges_to_add:
        if e[0] not in grid:
            print(f"{e[0]} is not in the grid")
        if e[1] not in grid:
            print(f"{e[1]} is not in the grid")
        grid.add_edge(e[0], e[1], weight=cost_shape_edge, crossing=[])

    # Similar to above, but easier, since we go over the shape and adapting the grid, so no tracking needed.
    for u, v in shape.edges():
        if (u, v) not in crossings:
            # Retain uncrossed shape edges, which are not yet in the graph
            grid.add_edge(u, v, weight=cost_shape_edge, crossing=[])
            continue
        edge_crossings = crossings[(u, v)]
        edge_crossings.sort(key=lambda x: dist(IDs[x][0], IDs[x][1], shape.nodes()[u]["x"], shape.nodes()[u]["y"]))
        for c in edge_crossings:
            if c not in grid:
                grid.add_node(c, x=IDs[c][0], y=IDs[c][1])
        for i in range(len(edge_crossings) - 1):
            grid.add_edge(edge_crossings[i], edge_crossings[i + 1], weight=cost_shape_edge, crossing=[])
        grid.add_edge(u, edge_crossings[0], weight=cost_shape_edge, crossing=[])
        grid.add_edge(edge_crossings[len(edge_crossings) - 1], v, weight=cost_shape_edge, crossing=[])

    # We mark nodes not on the shape as such
    for node, data in grid.nodes(data=True):
        if "shape_part" not in data or not data["shape_part"]:
            grid.nodes()[node]["shape_part"] = False
    return grid


# Removes (by contraction) the original vertices of the shape (which are of degree 2 and therefore not wanted or needed)
def remove_original_shape_nodes(shape, grid):
    to_contract = []
    for node in shape.nodes():
        if grid.degree[node] == 2:
            for n in grid.neighbors(node):
                nghb = n
                break
            to_contract.append((nghb, node))
    for n1, n2 in to_contract:
        if grid.has_node(n2) and grid.has_node(n1):
            # print(f"\t\t\tDEBUG: merging {n2} into {n1} because of DEGREE")
            nx.contracted_nodes(grid, n1, n2, self_loops=False, copy=False)
            # contraction is saved in the data dict and causes problems later on
            grid.nodes()[n1].pop("contraction")

    return grid


# We identify the cell in which a point lies, then based on the relative placement inside the cell, we add new nodes and edges NEW METHOD
def octolinearize_single_point_2(point, g_left, g_bot, g_width, g_height, cell_size, grid, max_ind, node):
    # double precision euqality check
    prec_bound = 0.000001

    x, y = point[0], point[1]
    x -= g_left
    y -= g_bot
    col = int((x + prec_bound) / cell_size)
    row = int((y + prec_bound) / cell_size)

    in_x = x % cell_size
    in_y = y % cell_size

    if in_x > cell_size - prec_bound:
        in_x -= cell_size

    if in_y > cell_size - prec_bound:
        in_y -= cell_size

    # input()

    # the index of the node, which is to the lower left of our node.
    i = col * g_height + row
    # print(f"called with height: {g_height} and width {g_width} resulting in col: {col} and row: {row} and therefore i: {i}")

    # cell offsets
    x_o, y_o = (col * cell_size) + g_left, (row * cell_size) + g_bot

    # In any case, we introduce 6 new points
    nodes_to_add = {}
    edges_to_add = []
    edges_to_remove = []

    edges_to_add_dict = {}

    if in_y > cell_size - prec_bound or in_x > cell_size - prec_bound:
        print(f"\t\t\t\t\t\t\t{node} exceeds expecations")

    W = cell_size
    v_h0 = (in_x, 0)
    v_lb = (in_x, in_x)
    v_rb = (in_x, W - in_x)
    v_h1 = (in_x, W)

    h_v0 = (0, in_y)
    h_lb = (in_y, in_y)
    h_rb = (W - in_y, in_y)
    h_v1 = (W, in_y)

    lb_v0 = (0, in_y - in_x)
    lb_h0 = (in_x - in_y, 0)
    lb_rb = ((in_x - in_y + W) / 2, (in_y - in_x + W) / 2)
    lb_v1 = (W, W - (in_x - in_y))
    lb_h1 = (W - (in_y - in_x), W)

    rb_v0 = (0, in_x + in_y)
    rb_h0 = (in_x + in_y, 0)
    rb_lb = ((in_x + in_y) / 2, (in_x + in_y) / 2)
    rb_v1 = (W, in_x + in_y - W)
    rb_h1 = (in_x + in_y - W, W)

    # Offsets
    offsets = (x_o, y_o)
    v_h0 = tuple(sum(a) for a in zip(v_h0, offsets))
    v_lb = tuple(sum(a) for a in zip(v_lb, offsets))
    v_rb = tuple(sum(a) for a in zip(v_rb, offsets))
    v_h1 = tuple(sum(a) for a in zip(v_h1, offsets))
    h_v0 = tuple(sum(a) for a in zip(h_v0, offsets))
    h_lb = tuple(sum(a) for a in zip(h_lb, offsets))
    h_rb = tuple(sum(a) for a in zip(h_rb, offsets))
    h_v1 = tuple(sum(a) for a in zip(h_v1, offsets))
    lb_v0 = tuple(sum(a) for a in zip(lb_v0, offsets))
    lb_h0 = tuple(sum(a) for a in zip(lb_h0, offsets))
    lb_rb = tuple(sum(a) for a in zip(lb_rb, offsets))
    lb_v1 = tuple(sum(a) for a in zip(lb_v1, offsets))
    lb_h1 = tuple(sum(a) for a in zip(lb_h1, offsets))
    rb_v0 = tuple(sum(a) for a in zip(rb_v0, offsets))
    rb_h0 = tuple(sum(a) for a in zip(rb_h0, offsets))
    rb_lb = tuple(sum(a) for a in zip(rb_lb, offsets))
    rb_v1 = tuple(sum(a) for a in zip(rb_v1, offsets))
    rb_h1 = tuple(sum(a) for a in zip(rb_h1, offsets))

    H = g_height

    l_v0 = (i, i + 1)
    l_v1 = (i + H, i + H + 1)
    l_h0 = (i, i + H)
    l_h1 = (i + 1, i + H + 1)
    l_lb = (i, i + H + 1)
    l_rb = (i + 1, i + H)

    if in_x > in_y:
        # Below bl-tr diagonal
        # print(f"{node} is below bl-tr diagonal")
        if in_x + in_y < W:
            # Below br-tl diagonal
            # print(f"{node} is below br-tl diagonal")
            #  /\
            # /  \
            # ----
            p1 = v_h0
            p2 = rb_h0
            p3 = h_rb
            p4 = lb_rb
            if in_x < W / 2:
                p5 = v_lb
            else:
                p5 = v_rb
            p6 = rb_lb
            p7 = h_lb
            p8 = lb_h0

            nodes_to_add[max_ind + 1] = p1
            nodes_to_add[max_ind + 2] = p2
            nodes_to_add[max_ind + 3] = p3
            nodes_to_add[max_ind + 4] = p4
            nodes_to_add[max_ind + 5] = p5
            nodes_to_add[max_ind + 6] = p6
            nodes_to_add[max_ind + 7] = p7
            nodes_to_add[max_ind + 8] = p8

            edges_to_remove.append(l_h0)
            # always have to be removed
            edges_to_remove.append(l_rb)
            edges_to_remove.append(l_lb)

            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_rb, max_ind + 3)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_rb, max_ind + 4)
            if in_x < W / 2:
                remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_lb, max_ind + 5)
            else:
                remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_rb, max_ind + 5)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_lb, max_ind + 6)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_lb, max_ind + 7)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_h0, max_ind + 8)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_h0, max_ind + 1)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_h0, max_ind + 2)
        else:
            # Above br-tl diagonal
            # print(f"{node} is above br-tl diagonal")
            #  /|
            #  \|
            p1 = v_rb
            p2 = rb_v1
            p3 = h_v1
            p4 = lb_v1
            p5 = v_lb
            p6 = rb_lb
            if in_y < W / 2:
                p7 = h_rb
            else:
                p7 = h_lb
            p8 = lb_rb

            nodes_to_add[max_ind + 1] = p1
            nodes_to_add[max_ind + 2] = p2
            nodes_to_add[max_ind + 3] = p3
            nodes_to_add[max_ind + 4] = p4
            nodes_to_add[max_ind + 5] = p5
            nodes_to_add[max_ind + 6] = p6
            nodes_to_add[max_ind + 7] = p7
            nodes_to_add[max_ind + 8] = p8

            edges_to_remove.append(l_v1)
            # always have to be removed
            edges_to_remove.append(l_rb)
            edges_to_remove.append(l_lb)

            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_rb, max_ind + 8)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_rb, max_ind + 1)
            if in_y < W / 2:
                remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_rb, max_ind + 7)
            else:
                remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_lb, max_ind + 7)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_lb, max_ind + 5)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_lb, max_ind + 6)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_v1, max_ind + 2)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_v1, max_ind + 3)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_v1, max_ind + 4)
    else:
        # Above bl-tr diagonal
        # print(f"{node} is above bl-tr diagonal")
        if in_x + in_y > W:
            # Above br-tl diagonal
            # print(f"{node} is above br-tl diagonal")
            # ----
            # \  /
            #  \/
            if in_x > W / 2:
                p1 = v_lb
            else:
                p1 = v_rb
            p2 = rb_lb
            p3 = h_lb
            p4 = lb_h1
            p5 = v_h1
            p6 = rb_h1
            p7 = h_rb
            p8 = lb_rb

            nodes_to_add[max_ind + 1] = p1
            nodes_to_add[max_ind + 2] = p2
            nodes_to_add[max_ind + 3] = p3
            nodes_to_add[max_ind + 4] = p4
            nodes_to_add[max_ind + 5] = p5
            nodes_to_add[max_ind + 6] = p6
            nodes_to_add[max_ind + 7] = p7
            nodes_to_add[max_ind + 8] = p8

            edges_to_remove.append(l_h1)
            # always have to be removed
            edges_to_remove.append(l_rb)
            edges_to_remove.append(l_lb)

            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_rb, max_ind + 7)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_rb, max_ind + 8)
            if in_x > W / 2:
                remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_lb, max_ind + 1)
            else:
                remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_rb, max_ind + 1)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_lb, max_ind + 2)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_lb, max_ind + 3)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_h1, max_ind + 4)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_h1, max_ind + 5)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_h1, max_ind + 6)
        else:
            # Below br-tl diagonal
            # print(f"{node} is below br-tl diagonal")
            #  |\
            #  |/
            p1 = v_lb
            p2 = rb_lb
            if in_y < W / 2:
                p3 = h_lb
            else:
                p3 = h_rb
            p4 = lb_rb
            p5 = v_rb
            p6 = rb_v0
            p7 = h_v0
            p8 = lb_v0

            nodes_to_add[max_ind + 1] = p1
            nodes_to_add[max_ind + 2] = p2
            nodes_to_add[max_ind + 3] = p3
            nodes_to_add[max_ind + 4] = p4
            nodes_to_add[max_ind + 5] = p5
            nodes_to_add[max_ind + 6] = p6
            nodes_to_add[max_ind + 7] = p7
            nodes_to_add[max_ind + 8] = p8

            edges_to_remove.append(l_v0)
            # always have to be removed
            edges_to_remove.append(l_rb)
            edges_to_remove.append(l_lb)

            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_rb, max_ind + 4)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_rb, max_ind + 5)
            if in_y < W / 2:
                remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_lb, max_ind + 3)
            else:
                remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_rb, max_ind + 3)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_lb, max_ind + 1)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_lb, max_ind + 2)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_v0, max_ind + 6)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_v0, max_ind + 7)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, *l_v0, max_ind + 8)

    # END OF NEW STUFF

    # print(f"col: {col} and row: {row} and therefore i: {i}")
    # print(f"x was: {x} and y was: {y} resulting with cell size {cell_size} in in_x: {in_x} and in_y: {in_y}")

    for i in range(1, 9):
        edges_to_add.append((node, max_ind + i))
    max_ind += 8

    # for edge in edges_to_add:
    #     if edge[0] not in nodes_to_add:
    #         print(f"failed at {edge[0]} of {edge}")
    #         print(f"nodes: {nodes_to_add}")
    #         print()
    #         print(f"edges: {edges_to_add}")
    #         input()
    #     if edge[1] not in nodes_to_add:
    #         print(f"failed at {edge[1]} of {edge}")
    #         print(f"nodes: {nodes_to_add}")
    #         print()
    #         print(f"edges: {edges_to_add}")
    #         input()
    return nodes_to_add, edges_to_add, edges_to_remove, max_ind + 9, edges_to_add_dict


# Shoots octolinear rays from a vertex on the shape (the new ones presumably, but function uses the "shape_part" attribute) to the next edge.
# This should create too close nodes. These will be marked as not to be used for candidate sets (again via attribute, "candidate" = False)
def octolinearize_shape_nodes(grid, g_left, g_bot, g_width, g_height, cell_size, cost_shape_connection_edge):
    nodes_to_add = {}
    edges_to_add_dict = {}
    edges_to_add = []
    edges_to_remove = []
    # All introduced new points have to have a new index, so not to overwrite an existing node
    max_ind = 0
    for node in grid:
        if node > max_ind:
            max_ind = node

    for node, data in grid.nodes(data=True):
        if data["shape_part"]:
            # this node is part of the shape. Therefore it is at this point an intersection point with a singular edge of the grid
            # We are now octolinarizing this bit
            n2add, e2add, e2rem, max_ind, new_edge_dict = octolinearize_single_point_2((data['x'], data['y']), g_left,
                                                                                       g_bot, g_width, g_height,
                                                                                       cell_size, grid, max_ind, node)
            # print(f"octolinearizing {node}.\nFor adding {e2add}\nwe have to remove {e2rem}")
            # print(f"instead we use {new_edge_dict}")
            # joining two dicts
            nodes_to_add = {**nodes_to_add, **n2add}

            edges_to_add = edges_to_add + e2add
            edges_to_add_dict = {**edges_to_add_dict, **new_edge_dict}

            edges_to_remove = edges_to_remove + e2rem

    for node, pos in nodes_to_add.items():
        # print(f"adding: {node} at {pos}")
        grid.add_node(node, x=pos[0], y=pos[1], shape_part=False, not_candidate=True)
    for edge in edges_to_add:
        grid.add_edge(edge[0], edge[1], crossing=[], weight=cost_shape_connection_edge)
    for edge, ind_dist_pairs in edges_to_add_dict.items():
        # Sort entries by distance to left endpoint
        ind_dist_pairs.sort(key=lambda x: x[1])
        grid.add_edge(edge[0], ind_dist_pairs[0][0], crossing=[], weight=cost_shape_connection_edge)
        for i in range(1, len(ind_dist_pairs)):
            # print(f"DEBUG: Adding {(ind_dist_pairs[i-1][0], ind_dist_pairs[i][0])}")
            grid.add_edge(ind_dist_pairs[i - 1][0], ind_dist_pairs[i][0], crossing=[],
                          weight=cost_shape_connection_edge)
        grid.add_edge(ind_dist_pairs[-1][0], edge[1], crossing=[], weight=cost_shape_connection_edge)
    for edge in edges_to_remove:
        if grid.has_edge(edge[0], edge[1]):
            grid.remove_edge(edge[0], edge[1])
