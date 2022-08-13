import math

from ac_metro.utility.abstractUtility import AbstractUtility
from ac_metro.utility.geometry.boundingBoxUtility import BoundingBoxUtility


class GridRasterizer(AbstractUtility):

    def execute(map, options=None):

        raster_resolution = options["resolution"] if "resolution" in options else 1

        # Compute Bounding Box properties
        if "bb_width" not in options or "bb_height" not in options:
            min_x, min_y, bb_width, bb_height = BoundingBoxUtility.execute(map, {"padding": options["padding"] if "padding" in options else 0})
        else:
            bb_width, bb_height = options["bb_width"], options["bb_height"]

        bb_area = bb_width * bb_height
        bb_aspect_ratio = bb_width / bb_height

        # X * Y = temp = ⌈A / (f * D)²⌉ from the paper
        temp = math.ceil(bb_area / math.pow(raster_resolution, 2))

        # substituting X in X*Y = temp -> Y²*r = temp -> Y = √(temp/r)
        Y = math.sqrt(temp / bb_aspect_ratio)

        # X, Y are still in the same ratio as the BB of the input
        X = Y * bb_aspect_ratio

        # grid size is integer
        grid_rows, grid_columns = math.ceil(X) + 1, math.ceil(Y) + 1

        return grid_rows, grid_columns
