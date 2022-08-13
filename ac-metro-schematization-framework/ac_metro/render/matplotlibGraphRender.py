import matplotlib.pyplot as plt
from matplotlib import collections  as mc

from ac_metro.render.abstractRender import AbstractRender

class MatplotlibGraphRender(AbstractRender):
    '''
    Use matplotlib to render a metro map as a simple graph.
    '''

    def render(map, options=None):
        fig, ax = plt.subplots()
        
        segments = []
        for u,v in map.graph.edges():
            
            st1 = map.graph.nodes[u]['station']
            st2 = map.graph.nodes[v]['station']

            segments.append([(st1.x, st1.y), (st2.x, st2.y)])

        lc = mc.LineCollection(segments, colors='#000000')
        ax.add_collection(lc)

        x = []
        y = []
        for v,data in map.graph.nodes(data=True):
            station = data['station']
            x.append(station.x)
            y.append(station.y)

        ax.scatter(x,y)

        ax.axes.set_aspect('equal')
        plt.show()