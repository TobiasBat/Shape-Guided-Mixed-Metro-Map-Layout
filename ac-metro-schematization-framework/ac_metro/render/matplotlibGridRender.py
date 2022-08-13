import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib import collections as mc
import os

from ac_metro.render.abstractRender import AbstractRender


class MatplotlibGridRender(AbstractRender):
    '''
    Use matplotlib to render the map.
    '''

    @staticmethod
    def render(map, options=None):
        fig_size = 15
        if "fig_size" in options:
            fig_size = options["fig_size"]
        fig, ax = plt.subplots(figsize=(fig_size, fig_size))

        dpi_value = 150
        if "dpi" in options:
            dpi_value = options["dpi"]

        grid = options["grid"]
        edge_colors = ["grey"] * len(grid.edges())
        pos = {}
        for node, data in grid.nodes(data=True):
            # print(f"DEBUG: {node} has {data}")
            pos[node] = np.array([data["x"], data["y"]])

        # for i, (u, v, data) in enumerate(grid.edges(data=True)):
        #     if math.isinf(data['weight']):
        #         edge_colors[i] = "deeppink"

        # for i, (u, v, data) in enumerate(grid.edges(data=True)):
        #     if grid.nodes()[u]["shape_part"] and grid.nodes()[v]["shape_part"]:
        #         edge_colors[i] = "deeppink"

        straight_prop = ""
        if "draw_straight" in options:
            straight_prop = options["draw_straight"]
            print(f"drawing {straight_prop} straight")

        if "highlight_nodes" in options:
            for i in options["highlight_nodes"]:
                center = (grid.nodes()[i]["x"], grid.nodes()[i]["y"])
                circle = plt.Circle(center, 2, color='deeppink')
                plt.gca().add_patch(circle)

        if "font_size" in options:
            font_size = options["font_size"]
        else:
            font_size = 10

        grid_labels = "grid_labels" in options and options["grid_labels"]

        if "show_nodes" in options and options["show_nodes"]:
            node_size = 1
        else:
            node_size = 0

        alpha_value = 0.2
        width_value = 1
        if "just_grid" in options and options["just_grid"]:
            alpha_value = 1
            width_value = 0.5
            node_size = 0.2
        nx.draw_networkx(grid, pos, node_size=node_size, edge_color=edge_colors, alpha=alpha_value, width=width_value,
                         with_labels=grid_labels, font_size=font_size)
        if "grid_edge_labels" in options:
            if options["grid_edge_labels"] == "weight":
                edge_labels = {key: f"{key}:{round(value, ndigits=2)}" for key, value in
                               nx.get_edge_attributes(grid, 'weight').items()}
            elif options["grid_edge_labels"] == "crossing":
                edge_labels = {key: f"{key}:{value}" for key, value in
                               nx.get_edge_attributes(grid, 'crossing').items()}

            nx.draw_networkx_edge_labels(grid, pos, alpha=1, edge_labels=edge_labels, font_size=1)

        if "just_grid" not in options or not options["just_grid"]:
            for key, stations in map.getLines():

                segments = []

                if map.graph.has_edge(stations[-1], stations[0]):
                    last = stations[-1]
                    start_at = 0
                else:
                    last = stations[0]
                    start_at = 1
                for s in stations[start_at:]:
                    st1 = map.graph.nodes[last]['station']
                    st2 = map.graph.nodes[s]['station']

                    c_points = map.graph[last][s]["edge"].controlPoints
                    # print(f"this is edge data {map.graph[last][s]}")
                    if len(c_points) == 0 or (
                            straight_prop != "" and straight_prop in map.graph[last][s] and map.graph[last][s][
                        straight_prop]):
                        # print("just first to last")
                        segments.append([(st1.x, st1.y), (st2.x, st2.y)])
                    else:
                        if map.graph[last][s]["edge"].source == last:
                            # print("OI!")
                            segments.append([(st1.x, st1.y), (c_points[0][0], c_points[0][1])])
                            last_point = c_points[0]
                            for point in c_points[1:]:
                                segments.append([(last_point[0], last_point[1]), (point[0], point[1])])
                                last_point = point
                            segments.append([(last_point[0], last_point[1]), (st2.x, st2.y)])
                        else:
                            # print("This aint it chief")
                            segments.append([(st2.x, st2.y), (c_points[0][0], c_points[0][1])])
                            last_point = c_points[0]
                            for point in c_points[1:]:
                                segments.append([(last_point[0], last_point[1]), (point[0], point[1])])
                                last_point = point
                            segments.append([(last_point[0], last_point[1]), (st1.x, st1.y)])

                    last = s

                color = map.lineColors[key] if key in map.lineColors else "#000000"

                lc = mc.LineCollection(segments, colors=color)
                ax.add_collection(lc)

            x = []
            y = []
            id_s = []
            for v, data in map.graph.nodes(data=True):
                station = data['station']
                x.append(station.x)
                y.append(station.y)
                id_s.append(station.id)

            node_colors = []

            for node in map.graph.nodes():
                if "settled_nodes" in options:
                    if node in options["settled_nodes"]:
                        node_colors.append("green")
                    else:
                        node_colors.append("black")
                else:
                    node_colors.append("black")

            ax.scatter(x, y, c=node_colors)

            if "station_names" in options and options["station_names"]:
                for i, txt in enumerate(id_s):
                    ax.annotate(txt, (x[i], y[i]))

        ax.axes.set_aspect('equal')

        if "to_file" in options:
            head = os.path.split(options["to_file"])[0]
            if not os.path.exists(str(head)):
                os.makedirs(str(head))
            plt.savefig(options["to_file"], dpi=dpi_value)
            plt.close()
        else:
            plt.show()
