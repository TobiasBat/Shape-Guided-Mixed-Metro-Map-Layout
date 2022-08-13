from ac_metro.utility.abstractUtility import AbstractUtility
import networkx as nx


class GridCreator(AbstractUtility):

    def execute(map, options=None):
        g = None
        x_offset = options["bottom_left"][0] if "bottom_left" in options else 0
        y_offset = options["bottom_left"][1] if "bottom_left" in options else 0
        grid_cost = options["grid_cost"] if "grid_cost" in options else 1
        grid_diagonal_cost = options["grid_diagonal_cost"] if "grid_diagonal_cost" in options else grid_cost

        if options["type"] == "octolinear":
            g = GridCreator.create_rectilinear_grid(options["width"], options["height"], options["x_scale"],
                                                    options["y_scale"], diagonals=True, grid_cost=grid_cost,
                                                    grid_diagonal_cost=grid_diagonal_cost, x_off=x_offset, y_off=y_offset)
        elif options["type"] == "rectilinear":
            g = GridCreator.create_rectilinear_grid(options["width"], options["height"], options["x_scale"],
                                                    options["y_scale"], diagonals=False, grid_cost=grid_cost,
                                                    grid_diagonal_cost=grid_diagonal_cost)
        elif options["type"] == "triangular":
            g = GridCreator.create_regular_grid(options["width"], options["height"], options["x_scale"],
                                                options["y_scale"])
        return g

    @classmethod
    def create_rectilinear_grid(cls, width, height, x_scale, y_scale, diagonals, grid_cost, grid_diagonal_cost, x_off=0, y_off=0):
        print(f"creating grid with normal: {grid_cost} and diag: {grid_diagonal_cost}")

        n_id = 0
        g = nx.Graph()
        # crossings = set()

        # print(f"creating grid of height {height} and width {width}")
        for j in range(width):
            for i in range(height):
                # print(f"\tnode: {n_id}")
                g.add_node(n_id, x=(j * x_scale) + x_off, y=(i * y_scale) + y_off)

                if i > 0:
                    # print(f"\tedge: {n_id, n_id - 1}")
                    g.add_edge(n_id, n_id - 1, weight=grid_cost)
                if j > 0:
                    # print(f"\tedge: {n_id, n_id - height}")
                    g.add_edge(n_id, n_id - height, weight=grid_cost)
                if diagonals:
                    # crossings.add(((n_id, n_id - height - 1), (n_id, n_id - height + 1)))
                    if j > 0 and i > 0:
                        # print(f"\tedge: {n_id, n_id - height - 1}")
                        g.add_edge(n_id, n_id - height - 1, weight=grid_diagonal_cost)
                        # newlist = [(n_id, n_id - height + 1)]
                        # if newlist[0][0] == -1 or newlist[0][1] == -1:
                        #     print(newlist)
                        #     print(f"i: {i} and j: {j} bong")
                        #     input()
                        # g[n_id][n_id - height - 1]["crossing"] = newlist
                    if j > 0 and i < (height - 1):
                        # print(f"\tedge: {n_id, n_id - height + 1}")
                        g.add_edge(n_id, n_id - height + 1, weight=grid_diagonal_cost)
                        # newlist = [(n_id, n_id - height - 1)]
                        # if newlist[0][0] == -1 or newlist[0][1] == -1:
                        #     print(newlist)
                        #     print(f"i: {i} and j: {j} bing")
                        #     input()
                        # g[n_id][n_id - height + 1]["crossing"] = newlist
                n_id += 1
        return g

    @classmethod
    def create_regular_grid(cls, width, height, x_scale, y_scale):
        raise NotImplementedError
