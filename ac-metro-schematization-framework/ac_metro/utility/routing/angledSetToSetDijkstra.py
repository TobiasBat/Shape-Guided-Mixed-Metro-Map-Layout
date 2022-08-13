import math

from ac_metro.utility.abstractUtility import AbstractUtility
import numpy as np
import networkx as nx
from decimal import *


def dist(x1, y1, x2, y2):
    return math.sqrt(math.pow(x1-x2, 2)+math.pow(y1-y2, 2))


class AngledSetToSetDijkstra(AbstractUtility):
    '''
    Run Set to Set Dijkstra, which takes the angles between nodes into account.
    With different weights this can return a bend minimal path, distance minimal path or something in between.
    '''

    def execute(map, options=None):

        graph = options["graph"]
        occupied_nodes = options["occupied_nodes"]
        blocked_edges = options["occupied_edges"]
        src_pos, tgt_pos = options["src_pos"], options["tgt_pos"]
        debug = options["debug"]
        # print(blocked_edges)

        # We're adding two new nodes with guaranteed new labels to the graph, which will get removed at the ende again
        a = graph.number_of_nodes()
        while a in graph:
            a += 1
        b = a + 1
        while b in graph:
            b += 1
        srcs = options["srcs"]
        tgts = options["tgts"]
        graph.add_node(a)
        graph.add_node(b)
        for src in srcs:
            graph.add_edge(a, src, weight=0.01*dist(graph.nodes()[src]["x"], graph.nodes()[src]["y"], options["src_coords"][0], options["src_coords"][1]))
            # graph.add_edge(a, src, weight=1)
            # DEBUG(debug, f"added {a}-{src}")
            # DEBUG(debug, f"a: {graph.nodes[a]}")
            # DEBUG(debug, f"src: {graph.nodes[src]}")
        for tgt in tgts:
            graph.add_edge(tgt, b, weight=0.01*dist(graph.nodes()[tgt]["x"], graph.nodes()[tgt]["y"], options["tgt_coords"][0], options["tgt_coords"][1]))
            # graph.add_edge(tgt, b, weight=1)
            # DEBUG(debug, f"added {tgt}-{b}")
        factor = options["factor"]
        path = AngledSetToSetDijkstra.dijkstra_angle(graph, a, b, factor, occupied_nodes, blocked_edges, src_pos, tgt_pos, debug)
        graph.remove_node(a)
        graph.remove_node(b)
        # print(f"we removed {a} and {b} and the path is currently {path}")
        # input()
        path = path[1:-1]
        path.reverse()
        DEBUG(debug, f"after handling it looks like this {path}")
        return path

    @staticmethod
    def dijkstra_angle(G: nx.Graph, src, tgt, factor, occupied, blocked, src_pos, tgt_pos, debug):
        n = len(G.nodes())
        cost = [100000] * n
        predList = [-1] * n
        cost[src] = 0
        queue = [src]
        
        DEBUG(debug, f"src position: {src_pos}")
        DEBUG(debug, f"tgt position: {tgt_pos}")
        # input()

        # for n in G.neighbors(src):
        #     c = cost[src] + G[src][n]['weight']
        #     if (c < cost[n]):
        #         if predList[n] < 0:
        #             queue.append(n)
        #
        #         cost[n] = c
        #         predList[n] = src

        debug_step = debug

        while True:
            queue.sort(reverse=True, key=lambda x: cost[x])
            # print(f"{queue}")
            DEBUG(debug_step, f"\nqueue at start: {queue}")
            DEBUG(debug_step, f"with these costs:")
            for node in queue:
                DEBUG(debug_step, f"{node} cost: {cost[node]}")
            best = queue.pop()
            # print(f"{best} is in the graph: {G.has_node(best)}")
            # Here we could break if best is the target node for a source to target dijkstra
            if predList[best] == -1:
                pred = None
            else:
                pred = G.nodes[predList[best]]
            bestN = G.nodes[best]

            DEBUG(debug_step, f"best in queue: {best} at {cost[best]}")
            # DEBUG(debug, f"neighbors of {best} are:")
            # for n in G.neighbors(best):
            #     DEBUG(debug, f"n: {n}")

            for n in G.neighbors(best):

                next = G.nodes[n]
                
                # THESE CHECKS ARE NOW INCORPORATED INTO THE EDGE WEIGHTS OF THE GRID, WHICH ARE UPDATED AFTER EVERY DIJKSTRA CALL
                # if n in occupied and src_pos != n and tgt_pos != n:
                #     DEBUG(debug, f"Occupied: {occupied}")
                #     DEBUG(debug, f"src_pos: {src_pos}")
                #     DEBUG(debug, f"tgt_pos: {tgt_pos}")
                #     DEBUG(debug, "Broke, because of occupation")
                #     continue
                # elif (best, n) in blocked or (n, best) in blocked:
                #     DEBUG(debug, "Broke, because of blocked edge")
                #     continue
                # else:
                #     c = cost[best] + G[best][n]['weight'] + AngledSetToSetDijkstra.angleCost(pred, bestN, next, factor, predList[best], best, n)
                #

                c = cost[best] + G[best][n]['weight'] + AngledSetToSetDijkstra.angleCost(pred, bestN, next, factor, predList[best], best, n)

                # print(f"a neighbor of {best} is {n} and the cost is {c}. Previous cost was {cost[n]}")
                if c < cost[n]:
                    if predList[n] < 0:
                        # print(f"which is why we append it to the queue")
                        queue.append(n)

                    cost[n] = c
                    predList[n] = best
            if len(queue) <= 0:
                break

        path = [tgt]
        last = tgt


        while True:
            next = predList[last]
            path.append(next)

            if next == src:
                DEBUG(debug, f"Would return with {path}")
                if debug: input()
                return path
            if next < 0:
                DEBUG(debug, f"Would return with empty, while path is {path}")
                if debug: input()
                return []
            last = next

    @staticmethod
    def angleCost(pred, curr, next, factor, predID=-1, currID=-1, nextID=-1):
        if pred is None:
            # print(f"no predecessor exists for {curr}")
            return 0
        if 'x' not in pred:
            # print(f"predecessor {predID} is not embedded")
            return 0
        if 'x' not in curr:
            # print(f"current {currID} is not embedded")
            return 0
        if 'x' not in next:
            # print(f"successor {nextID} is not embedded")
            return 0
        predCoord = np.array([pred['x'], pred['y']])
        currCoord = np.array([curr['x'], curr['y']])
        nextCoord = np.array([next['x'], next['y']])

        uv = predCoord - currCoord
        vw = nextCoord - currCoord
        uvD = np.linalg.norm(uv)
        vwD = np.linalg.norm(vw)
        uv = uv / uvD
        vw = vw / vwD

        dot = round(np.dot(uv, vw), ndigits=6)
        ang = np.arccos(dot)

        if ang < 0:
            print('Somehow the cost for a bend was negative. This will wreak havoc with the Dijkstra results.')

        return factor * math.pow(np.pi - ang, 2)


def DEBUG(debug, msg):
    if debug:
        print(msg)