import matplotlib.pyplot as plt
from matplotlib import collections  as mc

from ac_metro.render.abstractRender import AbstractRender

def edge_to_segments(map, st1_id, st2_id):
    segments = []

    st1 = map.graph.nodes[st1_id]['station']
    st2 = map.graph.nodes[st2_id]['station']
    control_points = map.graph[st1_id][st2_id]["edge"].getControlPoints(st1_id)
    if len(control_points) == 0:
        segments.append([(st1.x, st1.y), (st2.x, st2.y)])
    else:
        segments.append([(st1.x, st1.y), (control_points[0][0], control_points[0][1])])
        last_point = control_points[0]
        for point in control_points[1:]:
            segments.append([(last_point[0], last_point[1]), (point[0], point[1])])
            last_point = point
        segments.append([(last_point[0], last_point[1]), (st2.x, st2.y)])

    return segments


class MatplotlibRender(AbstractRender):
    '''
    Use matplotlib to render the map.
    '''


    def render(map, options=None):
        fig, ax = plt.subplots()

        control_points = []
        show_control_points = False
        if 'show_control_points' in options:
            show_control_points = options['show_control_points']

        for key, stations in map.getLines():

            segments = []

            last = stations[0]
            for s in stations[1:]:
                segments.extend(edge_to_segments(map, last, s))
                control_points.extend(map.graph[last][s]["edge"].controlPoints)
                last = s

            if map.graph.has_edge(last, stations[0]):
                segments.extend(edge_to_segments(map, last, stations[0]))
                control_points.extend(map.graph[last][stations[0]]["edge"].controlPoints)

            lc = mc.LineCollection(segments, colors=(map.lineColors[key] if key in map.lineColors else '#000000'))
            ax.add_collection(lc)

        x = []
        y = []
        for v,data in map.graph.nodes(data=True):
            station = data['station']
            x.append(station.x)
            y.append(station.y)

        ax.scatter(x, y, c="blue")

        if show_control_points is True:
            control_points_x = []
            control_points_y = []
            for (control_point_x, control_point_y) in control_points:
                control_points_x.append(control_point_x)
                control_points_y.append(control_point_y)

            ax.scatter(control_points_x, control_points_y, c="red", marker="*")

        ax.axes.set_aspect('equal')
        plt.show()
