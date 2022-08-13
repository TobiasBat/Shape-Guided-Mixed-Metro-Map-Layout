from abc import ABC, abstractmethod

class AbstractExport(ABC):
    '''
    Abstract base class for exporting data
    '''

    @staticmethod
    @abstractmethod
    def export(options:dict):
        """
        This abstract method should export a map to a specific format. 
        
        :options: Optional dictionary with subclass specific options. 

        :return: None
        """
        raise NotImplementedError()