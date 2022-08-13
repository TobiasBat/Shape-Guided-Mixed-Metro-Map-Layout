import math

from ac_metro.utility.abstractUtility import AbstractUtility


class GridBoundingBox(AbstractUtility):

    # Computes the bounding box of the grid
    @staticmethod
    def execute(map, options=None):
        grid = options["grid"] if "grid" in options else print(f"No grid provided to {GridBoundingBox.__name__}")
        pad = options["padding"] if "padding" in options else 0
        min_x, min_y, max_x, max_y = math.inf, math.inf, -math.inf, -math.inf
        for i, (node, data) in enumerate(grid.nodes(data=True)):
            min_x = min(min_x, data["x"])
            max_x = max(max_x, data["x"])
            min_y = min(min_y, data["y"])
            max_y = max(max_y, data["y"])

        return min_x - pad, min_y - pad, (max_x - min_x) + (2 * pad), (max_y - min_y) + (2 * pad)