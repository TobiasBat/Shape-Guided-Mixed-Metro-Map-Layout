import abc

from ac_metro.map.map import Map

class AbstractSchematic(abc.ABC):
    '''
    Abstract base class to schematize a metro map. Assumes that each station has an x- and y-coordinate.
    '''

    @staticmethod
    @abc.abstractmethod
    def schematize(map : Map, options=None):
        '''
        Schematize a map.

        :map: Map object that is schematized.
        :options: Optional options dictionary

        :return: None
        '''
        raise NotImplementedError