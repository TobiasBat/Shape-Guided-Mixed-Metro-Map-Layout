from ac_metro.utility.abstractUtility import AbstractUtility
from math import sin, cos

def rotate_point(cx, cy, angle, sx, sy):
    s = sin(angle)
    c = cos(angle)
    sx -= cx
    sy -= cy

    xnew = sx * c - sy * s
    ynew = sx * s + sy * c

    sx = xnew + cx
    sy = ynew + cy
    return sx, sy

class RotateUtility(AbstractUtility):
    '''
    Use a dictionary with identifier-color pairs to color a metro map.
    '''

    def execute(map, options=None):
        if options == None:
            return

        if "rotate_centroid" not in options:
            print("No centroid provided to RotateUtility")
            return
        else:
            centroid = options["rotate_centroid"]
            c_x = centroid[0]
            c_y = centroid[1]

        if "rotate_angle" not in options:
            print("No angle provided to RotateUtility")
            return
        else:
            angle = options["rotate_angle"]

        for i, (node, data) in enumerate(map.graph.nodes(data=True)):
            station = data["station"]
            station.x, station.y = rotate_point(c_x, c_y, angle, station.x, station.y)