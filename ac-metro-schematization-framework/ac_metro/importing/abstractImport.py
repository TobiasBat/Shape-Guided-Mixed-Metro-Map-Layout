from abc import ABC, abstractmethod

class AbstractImport(ABC):
    '''
    Abstract base class for importing data to the frameworks used format
    '''

    @staticmethod
    @abstractmethod
    def load(options:dict):
        """
        This abstract method should return a map object. Options are specified as an dictionary. 
        
        :options: Optional dictionary with subclass specific options. 

        :return: Map
        """
        raise NotImplementedError()