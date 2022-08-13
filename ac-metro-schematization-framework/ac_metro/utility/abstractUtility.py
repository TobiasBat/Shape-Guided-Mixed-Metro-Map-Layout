import abc

class AbstractUtility(abc.ABC):
    '''
    Abstract base class to provide utility functionality.  (e.g. Line coloring, Scaling)
    '''

    @staticmethod
    @abc.abstractmethod
    def execute(map, options=None):
        '''
        Execute the utility function

        :map: Map object
        :options: Optional options dictionary

        :return: None
        '''
        raise NotImplementedError