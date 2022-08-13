from ac_metro.utility.abstractUtility import AbstractUtility
from ac_metro.utility.geometry.boundingBoxUtility import BoundingBoxUtility
from ac_metro.utility.grid.gridBoundingBoxUtility import GridBoundingBox
from ac_metro.utility.transform.mapTranslater import MapTranslater


class MapCenterer(AbstractUtility):

    @staticmethod
    def execute(map, options=None):
        grid = options["grid"]
        g_bb_left, g_bb_bot, g_bb_width, g_bb_height = GridBoundingBox.execute(None, {"grid": grid})
        bb_left, bb_bot, bb_width, bb_height = BoundingBoxUtility.execute(map, {"padding": 0})

        x_diff = (g_bb_width - bb_width) / 2
        y_diff = (g_bb_height - bb_height) / 2

        # Sanity Check
        if x_diff < 0 or y_diff < 0:
            print(f"The grid has a smaller bounding box than the graph!")

        t_vector = (g_bb_left - bb_left + x_diff, g_bb_bot - bb_bot + y_diff)
        MapTranslater.execute(map, {"translate_vector": t_vector})
