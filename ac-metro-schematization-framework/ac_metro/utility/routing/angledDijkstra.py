from ac_metro.utility.abstractUtility import AbstractUtility
import numpy as np
import networkx as nx

class AngledDijkstra(AbstractUtility):
    '''
    Run Dijkstra, which takes the angles between nodes into account.
    With different weights this can return a bend minimal path, distance minimal path or something in between.
    '''

    # Incorporate as provided argument

    def execute(map, options=None):
        graph = options["graph"]
        src = options["src"]
        tgt = options["tgt"]
        factor = options["factor"]
        return AngledDijkstra.dijkstra_angle(graph, src, tgt, factor)

    @staticmethod
    def dijkstra_angle(G: nx.Graph, src, tgt, factor):
        n = len(G.nodes())
        cost = [100000] * n
        predList = [-1] * n
        cost[src] = 0
        queue = []

        for n in G.neighbors(src):
            c = cost[src] + G[src][n]['weight']
            if (c < cost[n]):
                if predList[n] < 0:
                    queue.append(n)

                cost[n] = c
                predList[n] = src

        while len(queue) > 0:
            queue.sort(reverse=True, key=lambda x: cost[x])
            best = queue.pop()
            pred = G.nodes[predList[best]]
            bestN = G.nodes[best]

            for n in G.neighbors(best):
                next = G.nodes[n]
                c = cost[best] + G[best][n]['weight'] + AngledDijkstra.angleCost(pred, bestN, next, factor)

                if c < cost[n]:
                    if predList[n] < 0:
                        queue.append(n)

                    cost[n] = c
                    predList[n] = best

        path = [tgt]
        last = tgt

        while True:
            print(last)
            next = predList[last]
            path.append(next)

            if next == src:
                return path
            if next < 0:
                return []
            last = next

    @staticmethod
    def angleCost(pred, curr, next, factor):
        predCoord = np.array([pred['x'], pred['y']])
        currCoord = np.array([curr['x'], curr['y']])
        nextCoord = np.array([next['x'], next['y']])

        uv = predCoord - currCoord
        vw = nextCoord - currCoord
        uvD = np.linalg.norm(uv)
        vwD = np.linalg.norm(vw)
        uv = uv / uvD
        vw = vw / vwD

        dot = np.dot(uv, vw)
        ang = np.arccos(dot)

        if ang < 0:
            print('Somehow the cost for a bend was negative. This will wreak havoc with the Dijkstra results.')

        return factor * (np.pi - ang)
