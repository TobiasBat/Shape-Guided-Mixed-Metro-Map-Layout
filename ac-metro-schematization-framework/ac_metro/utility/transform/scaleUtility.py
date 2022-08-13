import sys
from ac_metro.error.metroErrors import OptionsError
from ac_metro.map import Map
from ac_metro.utility.abstractUtility import AbstractUtility


class ScaleUtility(AbstractUtility):
    '''
    Use a dictionary with identifier-color pairs to color a metro map.
    '''

    def scale(map : Map, options):
        
        minX = sys.maxsize
        maxX = -sys.maxsize
        minY = sys.maxsize
        maxY = -sys.maxsize

        for station in map.getStations():
            x = station.x
            y = station.y

            minX = min(minX, x)
            maxX = max(maxX, x)
            minY = min(minY, y)
            maxY = max(maxY, y)

        # TODO add control points

        width = options['width']
        height = options['height']

        widthOld =  (maxX - minX)
        heightOld = (maxY - minY)
        
        if options["aspect"]:

            sNew = width / height
            sOld = widthOld / heightOld

            if sNew > sOld:
                xFactor = widthOld * height / heightOld
                yFactor = height
            else:
                xFactor = width
                yFactor = heightOld * width / widthOld

            for station in map.getStations():
                station.x = (station.x - minX) / (maxX - minX) * xFactor
                station.y = (station.y - minY) / (maxY - minY) * yFactor
            # TODO add control points
        else:
            for station in map.getStations():
                station.x = (station.x - minX) / (maxX - minX) * width
                station.y = (station.y - minY) / (maxY - minY) * height

            # TODO add control points

    def execute(map : Map, options=None):
        if options == None:
            raise OptionsError("ScaleUtility: No options given.")

        if "aspect" not in options:
            raise OptionsError("ScaleUtility: No aspects given.")

        if "width" not in options:
            raise OptionsError("ScaleUtility: No width given.")

        if "height" not in options:
            raise OptionsError("ScaleUtility: No height given.")

        ScaleUtility.scale(map, options)

        for snapshot in map.snapshots:
            ScaleUtility.scale(snapshot.map, options)