

# Shape handling
# TODO: Unfinished
def remove_close_stations(grid, shape, threshold):
    to_remove = []
    for g_node, g_data in grid.nodes(data=True):
        g_pos = (g_data["x"], g_data["y"])
        for s_node, s_data in shape.nodes(data=True):
            s_pos = (s_data["x"], s_data["y"])
            distance = dist(g_pos[0], g_pos[1], s_pos[0], s_pos[1])
            if distance < threshold:
                to_remove.append(g_node)
                break
    for node in to_remove:
        grid.remove_node(node)
        # ADD EDGES BACK IN BASED ON INDICES I GUESS

    for s_node, s_data in shape.nodes(data=True):
        s_x, s_y = s_data["x"], s_data["y"]
        for g_node, g_data in grid.nodes(data=True):
            g_x, g_y = g_data["x"], g_data["y"]
            distance = dist(s_x, s_y, g_x, g_y)
            if distance < threshold:
                custom_grid.add_edge(s_node, g_node, crossing=[], shape=True)


# Replaces crossing edges of shape and grid wih subdivided edges and nodes on the crossings
def intersections_to_nodes(grid, crossings, IDs, shape):
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
        grid.add_edge(e[0], e[1], weight=COST_SHAPE_EDGE, crossing=[])

    # Similar to above, but easier, since we go over the shape and adapting the grid, so no tracking needed.
    for u, v in shape.edges():
        if (u, v) not in crossings:
            # Retain uncrossed shape edges, which are not yet in the graph
            grid.add_edge(u, v, weight=COST_SHAPE_EDGE, crossing=[])
            continue
        edge_crossings = crossings[(u, v)]
        edge_crossings.sort(key=lambda x: dist(IDs[x][0], IDs[x][1], shape.nodes()[u]["x"], shape.nodes()[u]["y"]))
        for c in edge_crossings:
            if c not in grid:
                grid.add_node(c, x=IDs[c][0], y=IDs[c][1])
        for i in range(len(edge_crossings) - 1):
            grid.add_edge(edge_crossings[i], edge_crossings[i + 1], weight=COST_SHAPE_EDGE, crossing=[])
        grid.add_edge(u, edge_crossings[0], weight=COST_SHAPE_EDGE, crossing=[])
        grid.add_edge(edge_crossings[len(edge_crossings) - 1], v, weight=COST_SHAPE_EDGE, crossing=[])

    # We mark nodes not on the shape as such
    for node, data in grid.nodes(data=True):
        if "shape_part" not in data or not data["shape_part"]:
            grid.nodes()[node]["shape_part"] = False
    return grid


# Merges points in the grid (shape or grid) into other points in the grid, if they are very close
def contract_close_points_into_shape(shape, grid, threshhold):
    to_contract = {}
    for s_node, s_data in grid.nodes(data=True):
        if "shape_part" not in s_data or not s_data["shape_part"]:
            continue
        to_contract[s_node] = []
        s_x, s_y = s_data["x"], s_data["y"]
        for g_node, g_data in grid.nodes(data=True):
            if g_node == s_node:
                continue
            g_x, g_y = g_data["x"], g_data["y"]
            if dist(s_x, s_y, g_x, g_y) <= threshhold:
                to_contract[s_node].append(g_node)

    for n1, data in grid.nodes(data=True):
        if "shape_part" not in data or not data["shape_part"]:
            continue
        for n2 in to_contract[n1]:
            if grid.has_node(n2) and grid.has_node(n1):
                # print(f"\t\t\tDEBUG merging {n2} into {n1} because of DISTANCE")
                # print(f"\t\t\tDEBUG data\n{n1}: {grid.nodes()[n1]}\n{n2}: {grid.nodes()[n2]}")
                grid = nx.contracted_nodes(grid, n1, n2, self_loops=False)
                # contraction is saved in the data dict and causes problems later on
                grid.nodes()[n1].pop("contraction")

    return grid


# We identify the cell in which a point lies, then based on the relative placement inside the cell, we add new nodes and edges
def octolinearize_single_point(point, g_left, g_bot, g_width, g_height, cell_size, grid, max_ind, node):
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

    # the index of the node, which is to the lower left of our node.
    i = col * g_height + row

    # cell offsets
    x_o, y_o = (col * cell_size) + g_left, (row * cell_size) + g_bot

    # In any case, we introduce 6 new points
    p_1, p_2, p_3, p_4, p_5, p_6 = None, None, None, None, None, None

    nodes_to_add = {}
    edges_to_add = []
    edges_to_remove = []

    edges_to_add_dict = {}

    if in_y > cell_size - prec_bound or in_x > cell_size - prec_bound:
        print(f"\t\t\t\t\t\t\t{node} exceeds expecations")

    # Here we determine the correct positions and create nodes and edges in the grid depending on it
    if in_x < prec_bound:
        if in_y < cell_size / 2:
            # print("DEBUG: vert_b")
            p_1 = (-(cell_size - in_y) / 2, (cell_size + in_y) / 2)
            p_2 = ((cell_size - in_y) / 2, (cell_size + in_y) / 2)
            p_3 = (in_y, in_y)
            p_4 = (in_y / 2, in_y / 2)
            p_5 = (-           in_y / 2, in_y / 2)
            p_6 = (-           in_y, in_y)

            nodes_to_add[max_ind + 1] = (x_o + p_1[0], y_o + p_1[1])
            nodes_to_add[max_ind + 2] = (x_o + p_2[0], y_o + p_2[1])
            nodes_to_add[max_ind + 3] = (x_o + p_3[0], y_o + p_3[1])
            nodes_to_add[max_ind + 4] = (x_o + p_4[0], y_o + p_4[1])
            nodes_to_add[max_ind + 5] = (x_o + p_5[0], y_o + p_5[1])
            nodes_to_add[max_ind + 6] = (x_o + p_6[0], y_o + p_6[1])

            edges_to_add.append((i, max_ind + 3))
            edges_to_add.append((max_ind + 3, max_ind + 4))
            edges_to_add.append((max_ind + 4, i + 1 + g_height))
            edges_to_remove.append((i, i + 1 + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + 1 + g_height, max_ind + 3)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + 1 + g_height, max_ind + 4)

            edges_to_add.append((i + 1, max_ind + 2))
            edges_to_add.append((max_ind + 2, i + g_height))
            edges_to_remove.append((i + 1, i + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + 1, i + g_height, max_ind + 2)

            edges_to_add.append((i, max_ind + 6))
            edges_to_add.append((max_ind + 6, max_ind + 5))
            edges_to_add.append((max_ind + 5, i - g_height + 1))
            edges_to_remove.append((i, i - g_height + 1))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i - g_height + 1, i, max_ind + 5)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i - g_height + 1, i, max_ind + 6)

            edges_to_add.append((i - g_height, max_ind + 1))
            edges_to_add.append((max_ind + 1, i + 1))
            edges_to_remove.append((i - g_height, i + 1))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i - g_height, i + 1, max_ind + 1)

        else:
            # print("DEBUG: vert_t")
            p_1 = (-(cell_size - in_y) / 2, (in_y + cell_size) / 2)
            p_2 = ((cell_size - in_y) / 2, (in_y + cell_size) / 2)
            p_3 = (cell_size - in_y, in_y)
            p_4 = (in_y / 2, in_y / 2)
            p_5 = (-           in_y / 2, in_y / 2)
            p_6 = (-(cell_size - in_y), in_y)

            nodes_to_add[max_ind + 1] = (x_o + p_1[0], y_o + p_1[1])
            nodes_to_add[max_ind + 2] = (x_o + p_2[0], y_o + p_2[1])
            nodes_to_add[max_ind + 3] = (x_o + p_3[0], y_o + p_3[1])
            nodes_to_add[max_ind + 4] = (x_o + p_4[0], y_o + p_4[1])
            nodes_to_add[max_ind + 5] = (x_o + p_5[0], y_o + p_5[1])
            nodes_to_add[max_ind + 6] = (x_o + p_6[0], y_o + p_6[1])

            edges_to_add.append((i + 1, max_ind + 2))
            edges_to_add.append((max_ind + 2, max_ind + 3))
            edges_to_add.append((max_ind + 3, i + g_height))
            edges_to_remove.append((i + 1, i + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + 1, i + g_height, max_ind + 2)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + 1, i + g_height, max_ind + 3)

            edges_to_add.append((i, max_ind + 4))
            edges_to_add.append((max_ind + 4, i + 1 + g_height))
            edges_to_remove.append((i, i + 1 + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + 1 + g_height, max_ind + 4)

            edges_to_add.append((i - g_height, max_ind + 6))
            edges_to_add.append((max_ind + 6, max_ind + 1))
            edges_to_add.append((max_ind + 1, i + 1))
            edges_to_remove.append((i - g_height, i + 1))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i - g_height, i + 1, max_ind + 1)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i - g_height, i + 1, max_ind + 6)

            edges_to_add.append((i - g_height + 1, max_ind + 5))
            edges_to_add.append((max_ind + 5, i))
            edges_to_remove.append((i - g_height + 1, i))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i - g_height + 1, i + 1, max_ind + 5)

    elif in_y < prec_bound:
        if in_x < cell_size / 2:
            # print("DEBUG: hor_l")
            p_1 = (in_x / 2, in_x / 2)
            p_2 = (in_x, in_x)
            p_3 = ((in_x + cell_size) / 2, (cell_size - in_x) / 2)
            p_4 = ((in_x + cell_size) / 2, -(cell_size - in_x) / 2)
            p_5 = (in_x, -           in_x)
            p_6 = (in_x / 2, -           in_x / 2)

            nodes_to_add[max_ind + 1] = (x_o + p_1[0], y_o + p_1[1])
            nodes_to_add[max_ind + 2] = (x_o + p_2[0], y_o + p_2[1])
            nodes_to_add[max_ind + 3] = (x_o + p_3[0], y_o + p_3[1])
            nodes_to_add[max_ind + 4] = (x_o + p_4[0], y_o + p_4[1])
            nodes_to_add[max_ind + 5] = (x_o + p_5[0], y_o + p_5[1])
            nodes_to_add[max_ind + 6] = (x_o + p_6[0], y_o + p_6[1])

            edges_to_add.append((i, max_ind + 1))
            edges_to_add.append((max_ind + 1, max_ind + 2))
            edges_to_add.append((max_ind + 2, i + 1 + g_height))
            edges_to_remove.append((i, i + 1 + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + 1 + g_height, max_ind + 1)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + 1 + g_height, max_ind + 2)

            edges_to_add.append((i + 1, max_ind + 3))
            edges_to_add.append((max_ind + 3, i + g_height))
            edges_to_remove.append((i + 1, i + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + 1, i + g_height, max_ind + 3)

            edges_to_add.append((i, max_ind + 6))
            edges_to_add.append((max_ind + 6, max_ind + 5))
            edges_to_add.append((max_ind + 5, i - 1 + g_height))
            edges_to_remove.append((i, i - 1 + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i - 1 + g_height, max_ind + 5)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i - 1 + g_height, max_ind + 6)

            edges_to_add.append((i + g_height, max_ind + 4))
            edges_to_add.append((max_ind + 4, i - 1))
            edges_to_remove.append((i + g_height, i - 1))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i - 1, i + g_height, max_ind + 4)

        else:
            # print("DEBUG: hor_r")
            p_1 = (in_x / 2, in_x / 2)
            p_2 = (in_x, cell_size - in_x)
            p_3 = ((in_x + cell_size) / 2, (cell_size - in_x) / 2)
            p_4 = ((in_x + cell_size) / 2, -(cell_size - in_x) / 2)
            p_5 = (in_x, -(cell_size - in_x))
            p_6 = (in_x / 2, -           in_x / 2)

            nodes_to_add[max_ind + 1] = (x_o + p_1[0], y_o + p_1[1])
            nodes_to_add[max_ind + 2] = (x_o + p_2[0], y_o + p_2[1])
            nodes_to_add[max_ind + 3] = (x_o + p_3[0], y_o + p_3[1])
            nodes_to_add[max_ind + 4] = (x_o + p_4[0], y_o + p_4[1])
            nodes_to_add[max_ind + 5] = (x_o + p_5[0], y_o + p_5[1])
            nodes_to_add[max_ind + 6] = (x_o + p_6[0], y_o + p_6[1])

            edges_to_add.append((i + 1, max_ind + 2))
            edges_to_add.append((max_ind + 2, max_ind + 3))
            edges_to_add.append((max_ind + 3, i + g_height))
            edges_to_remove.append((i + 1, i + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + 1, i + g_height, max_ind + 2)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + 1, i + g_height, max_ind + 3)

            edges_to_add.append((i, max_ind + 1))
            edges_to_add.append((max_ind + 1, i + 1 + g_height))
            edges_to_remove.append((i, i + 1 + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + 1 + g_height, max_ind + 1)

            edges_to_add.append((i - 1, max_ind + 5))
            edges_to_add.append((max_ind + 5, max_ind + 4))
            edges_to_add.append((max_ind + 4, i + g_height))
            edges_to_remove.append((i - 1, i + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i - 1, i + g_height, max_ind + 4)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i - 1, i + g_height, max_ind + 5)

            edges_to_add.append((i, max_ind + 6))
            edges_to_add.append((max_ind + 6, i + g_height - 1))
            edges_to_remove.append((i, i + g_height - 1))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + g_height - 1, max_ind + 6)

    elif abs(in_x - in_y) < prec_bound:
        if in_x > cell_size / 2:
            # print("DEBUG: right_up_t")
            p_1 = (cell_size - in_x, in_y)
            p_2 = (-cell_size + 2 * in_x, cell_size)
            p_3 = (in_x, cell_size)
            p_4 = (cell_size, in_y)
            p_5 = (cell_size, -cell_size + 2 * in_y)
            p_6 = (in_x, cell_size - in_y)

            nodes_to_add[max_ind + 1] = (x_o + p_1[0], y_o + p_1[1])
            nodes_to_add[max_ind + 2] = (x_o + p_2[0], y_o + p_2[1])
            nodes_to_add[max_ind + 3] = (x_o + p_3[0], y_o + p_3[1])
            nodes_to_add[max_ind + 4] = (x_o + p_4[0], y_o + p_4[1])
            nodes_to_add[max_ind + 5] = (x_o + p_5[0], y_o + p_5[1])
            nodes_to_add[max_ind + 6] = (x_o + p_6[0], y_o + p_6[1])

            edges_to_add.append((i + 1, max_ind + 1))
            edges_to_add.append((max_ind + 1, max_ind + 6))
            edges_to_add.append((max_ind + 6, i + g_height))
            edges_to_remove.append((i + 1, i + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + 1, i + g_height, max_ind + 1)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + 1, i + g_height, max_ind + 6)

            edges_to_add.append((i + 1, max_ind + 2))
            edges_to_add.append((max_ind + 2, max_ind + 3))
            edges_to_add.append((max_ind + 3, i + 1 + g_height))
            edges_to_remove.append((i + 1, i + 1 + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + 1, i + 1 + g_height, max_ind + 2)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + 1, i + 1 + g_height, max_ind + 3)

            edges_to_add.append((i + g_height, max_ind + 5))
            edges_to_add.append((max_ind + 5, max_ind + 4))
            edges_to_add.append((max_ind + 4, i + 1 + g_height))
            edges_to_remove.append((i + g_height, i + 1 + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + g_height, i + 1 + g_height, max_ind + 4)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + g_height, i + 1 + g_height, max_ind + 5)

        else:
            # print("DEBUG: right_up_b")
            p_1 = (0, in_y)
            p_2 = (0, 2 * in_y)
            p_3 = (in_x, cell_size - in_y)
            p_4 = (cell_size - in_x, in_y)
            p_5 = (2 * in_x, 0)
            p_6 = (in_x, 0)

            nodes_to_add[max_ind + 1] = (x_o + p_1[0], y_o + p_1[1])
            nodes_to_add[max_ind + 2] = (x_o + p_2[0], y_o + p_2[1])
            nodes_to_add[max_ind + 3] = (x_o + p_3[0], y_o + p_3[1])
            nodes_to_add[max_ind + 4] = (x_o + p_4[0], y_o + p_4[1])
            nodes_to_add[max_ind + 5] = (x_o + p_5[0], y_o + p_5[1])
            nodes_to_add[max_ind + 6] = (x_o + p_6[0], y_o + p_6[1])

            edges_to_add.append((i, max_ind + 1))
            edges_to_add.append((max_ind + 1, max_ind + 2))
            edges_to_add.append((max_ind + 2, i + 1))
            edges_to_remove.append((i, i + 1))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + 1, max_ind + 1)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + 1, max_ind + 2)

            edges_to_add.append((i, max_ind + 6))
            edges_to_add.append((max_ind + 6, max_ind + 5))
            edges_to_add.append((max_ind + 5, i + g_height))
            edges_to_remove.append((i, i + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + g_height, max_ind + 5)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + g_height, max_ind + 6)

            edges_to_add.append((i + 1, max_ind + 3))
            edges_to_add.append((max_ind + 3, max_ind + 4))
            edges_to_add.append((max_ind + 4, i + g_height))
            edges_to_remove.append((i + 1, i + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + 1, i + g_height, max_ind + 3)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + 1, i + g_height, max_ind + 4)

    elif in_x + in_y > cell_size - prec_bound:
        if in_x < cell_size / 2:
            # print("DEBUG: left_up_t")
            p_1 = (0, in_y)
            p_2 = (in_x, cell_size)
            p_3 = (2 * in_x, cell_size)
            p_4 = (cell_size - in_x, in_y)
            p_5 = (in_x, cell_size - in_y)
            p_6 = (0, -cell_size + 2 * in_y)

            nodes_to_add[max_ind + 1] = (x_o + p_1[0], y_o + p_1[1])
            nodes_to_add[max_ind + 2] = (x_o + p_2[0], y_o + p_2[1])
            nodes_to_add[max_ind + 3] = (x_o + p_3[0], y_o + p_3[1])
            nodes_to_add[max_ind + 4] = (x_o + p_4[0], y_o + p_4[1])
            nodes_to_add[max_ind + 5] = (x_o + p_5[0], y_o + p_5[1])
            nodes_to_add[max_ind + 6] = (x_o + p_6[0], y_o + p_6[1])

            edges_to_add.append((i, max_ind + 6))
            edges_to_add.append((max_ind + 6, max_ind + 1))
            edges_to_add.append((max_ind + 1, i + 1))
            edges_to_remove.append((i, i + 1))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + 1, max_ind + 1)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + 1, max_ind + 6)

            edges_to_add.append((i + 1, max_ind + 2))
            edges_to_add.append((max_ind + 6, max_ind + 3))
            edges_to_add.append((max_ind + 3, i + 1 + g_height))
            edges_to_remove.append((i + 1, i + 1 + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + 1, i + 1 + g_height, max_ind + 3)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + 1, i + 1 + g_height, max_ind + 6)

            edges_to_add.append((i, max_ind + 5))
            edges_to_add.append((max_ind + 5, max_ind + 4))
            edges_to_add.append((max_ind + 4, i + 1 + g_height))
            edges_to_remove.append((i, i + 1 + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + 1 + g_height, max_ind + 4)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + 1 + g_height, max_ind + 5)

        else:
            # print("DEBUG: left_up_b")
            p_1 = (cell_size - in_x, in_y)
            p_2 = (in_x, cell_size - in_y)
            p_3 = (cell_size, 2 * in_y)
            p_4 = (cell_size, in_y)
            p_5 = (in_x, 0)
            p_6 = (-cell_size + 2 * in_x, 0)

            nodes_to_add[max_ind + 1] = (x_o + p_1[0], y_o + p_1[1])
            nodes_to_add[max_ind + 2] = (x_o + p_2[0], y_o + p_2[1])
            nodes_to_add[max_ind + 3] = (x_o + p_3[0], y_o + p_3[1])
            nodes_to_add[max_ind + 4] = (x_o + p_4[0], y_o + p_4[1])
            nodes_to_add[max_ind + 5] = (x_o + p_5[0], y_o + p_5[1])
            nodes_to_add[max_ind + 6] = (x_o + p_6[0], y_o + p_6[1])

            edges_to_add.append((i, max_ind + 1))
            edges_to_add.append((max_ind + 1, max_ind + 2))
            edges_to_add.append((max_ind + 2, i + 1 + g_height))
            edges_to_remove.append((i, i + 1 + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + 1 + g_height, max_ind + 1)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + 1 + g_height, max_ind + 2)

            edges_to_add.append((i, max_ind + 6))
            edges_to_add.append((max_ind + 6, max_ind + 5))
            edges_to_add.append((max_ind + 5, i + g_height))
            edges_to_remove.append((i, i + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + g_height, max_ind + 5)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i, i + g_height, max_ind + 6)

            edges_to_add.append((i + g_height, max_ind + 4))
            edges_to_add.append((max_ind + 4, max_ind + 3))
            edges_to_add.append((max_ind + 3, i + 1 + g_height))
            edges_to_remove.append((i + g_height, i + 1 + g_height))
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + g_height, i + 1 + g_height, max_ind + 3)
            remember_edge_to_add(grid, edges_to_add_dict, nodes_to_add, i + g_height, i + 1 + g_height, max_ind + 4)

    else:
        print("Crossing not on any line")

    # print(f"col: {col} and row: {row} and therefore i: {i}")
    # print(f"x was: {x} and y was: {y} resulting with cell size {cell_size} in in_x: {in_x} and in_y: {in_y}")

    for i in range(1, 7):
        edges_to_add.append((node, max_ind + i))
    max_ind += 6
    return nodes_to_add, edges_to_add, edges_to_remove, max_ind + 6, edges_to_add_dict
