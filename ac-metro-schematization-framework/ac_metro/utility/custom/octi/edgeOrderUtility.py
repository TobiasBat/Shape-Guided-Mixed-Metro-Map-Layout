from ac_metro.utility.abstractUtility import AbstractUtility
from random import shuffle, seed
import math

# Helper methods
def dist(x1, y1, x2, y2):
    '''
    Returns the Euclidean distance
    '''
    return math.sqrt(math.pow(x1-x2, 2)+math.pow(y1-y2, 2))

class EdgeOrderUtility(AbstractUtility):
    '''
    Compute an ordered set of all edges of the map based on a specified criteria
    '''

    @staticmethod
    def execute(map, options=None):
        if "seed" in options:
            seed(options["seed"])
            # print(f"seed set to {options['seed']}")

        method = options["method"]

        edge_list = []
        if "edge_list" in options:
            # print("We have an edge list")
            edge_list = options["edge_list"]
        else:
            for u, v in map.graph.edges():
                edge_list.append((u, v))

        lines_in_order = []

        if method == "line_degree":
            vertex_line_deg = {}
            for v in map.graph.nodes():
                value = 0
                for edge in map.graph.edges(v):
                    value += len(map.graph[edge[0]][edge[1]]["edge"].getLines())
                vertex_line_deg[v] = value
            for u, v in map.graph.edges():
                lines_in_order.append(((u, v), vertex_line_deg[u] + vertex_line_deg[v]))

        for u, v in edge_list:
            data = map.graph[u][v]

            # if method == "line_degree":
            #     value = len(data["edge"].getLines())
            if method == "length":
                s1, s2 = map.graph.nodes()[u]["station"], map.graph.nodes()[v]["station"]
                value = dist(s1.x, s1.y, s2.x, s2.y)
            elif method == "random":
                value = 0

            lines_in_order.append(((u, v), value))

        if method == "random":
            shuffle(lines_in_order)
        else:
            lines_in_order = sorted(lines_in_order, reverse=True, key=lambda x: x[1])
        
        return [i[0] for i in lines_in_order]

        #     for u, v in edge_list:
        #         data = map[u][v]
        #         lines_in_order.append(((u, v), len(data["edge"].getLines())))
        #     lines_in_order = sorted(lines_in_order, reverse=True, key=lambda x: x[1])
        #     return [i[0] for i in lines_in_order]

        # elif method == "length":
        #     for u, v in edge_list:
        #         data = map[u][v]
        #         s1, s2 = map.graph.nodes()[u]["station"], map.graph.nodes()[v]["station"]
        #         lines_in_order.append(((u, v), dist(s1.x, s1.y, s2.x, s2.y)))
        #     lines_in_order = sorted(lines_in_order, reverse=True, key=lambda x: x[1])
        #     return [i[0] for i in lines_in_order]

        # elif method == "random":
        #     for u, v, data in map.graph.edges(data=True):
        #         lines_in_order.append((u, v))
        #     shuffle(lines_in_order)
        #     return lines_in_order

        # if "prefer_att" in options:
        #     pref_att = options["prefer_att"]
        #     prefer_lines = options["preferred_edges"][0]
        #     rest_lines = options["preferred_edges"][1]

        #     if method == "line_degree":
        #         for u, v, data in map.graph.edges(data=True):
        #             if not map.graph[u][v]["edge"][pref_att]:
        #                 continue
        #             prefer_lines.append(((u, v), len(data["edge"].getLines())))
        #         lines_in_order = sorted(prefer_lines, reverse=True, key=lambda x: x[1])
        #         prefer_lines = [i[0] for i in lines_in_order]

        #         for u, v, data in map.graph.edges(data=True):
        #             if map.graph[u][v][pref_att]:
        #                 continue
        #             rest_lines.append(((u, v), len(data["edge"].getLines())))
        #         lines_in_order = sorted(rest_lines, reverse=True, key=lambda x: x[1])
        #         rest_lines = [i[0] for i in lines_in_order]

        #         prefer_lines.extend(rest_lines)
        #         return prefer_lines

        #     elif method == "length":
        #         for u, v, data in map.graph.edges(data=True):
        #             if not map.graph[u][v][pref_att]:
        #                 continue
        #             s1, s2 = map.graph.nodes()[u]["station"], map.graph.nodes()[v]["station"]
        #             prefer_lines.append(((u, v), dist(s1.x, s1.y, s2.x, s2.y)))
        #         lines_in_order = sorted(prefer_lines, reverse=True, key=lambda x: x[1])
        #         prefer_lines = [i[0] for i in lines_in_order]

        #         for u, v, data in map.graph.edges(data=True):
        #             if map.graph[u][v][pref_att]:
        #                 continue
        #             s1, s2 = map.graph.nodes()[u]["station"], map.graph.nodes()[v]["station"]
        #             rest_lines.append(((u, v), dist(s1.x, s1.y, s2.x, s2.y)))
        #         lines_in_order = sorted(rest_lines, reverse=True, key=lambda x: x[1])
        #         rest_lines = [i[0] for i in lines_in_order]

        #         prefer_lines.extend(rest_lines)
        #         return prefer_lines

        #     elif method == "random":
        #         for u, v, data in map.graph.edges(data=True):
        #             if not map.graph[u][v][pref_att]:
        #                 continue
        #             prefer_lines.append((u, v))
        #         shuffle(prefer_lines)

        #         for u, v, data in map.graph.edges(data=True):
        #             if map.graph[u][v][pref_att]:
        #                 continue
        #             rest_lines.append((u, v))
        #         shuffle(rest_lines)

        #         prefer_lines.extend(rest_lines)
        #         return prefer_lines

        # else:
