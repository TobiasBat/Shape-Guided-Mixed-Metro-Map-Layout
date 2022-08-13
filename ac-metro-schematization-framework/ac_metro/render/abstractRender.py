import abc

class AbstractRender(abc.ABC):
    '''
    Abstract class that specifies a common render interface.
    '''

    @staticmethod
    @abc.abstractmethod
    def render(map, options=None):
        '''
        Render a map. x- and y-coordinates are assumed to be present in the map.

        :map: The map object
        :options: Additional dictionary that specifies options that might be required for subclasses.

        :return: None
        '''
        raise NotImplementedError