from  ac_metro.importing.abstractImport import AbstractImport
import pickle

class PickleImport(AbstractImport):
    '''
    Load a metro map from an Graphml file.
    '''

    def load(options:dict):
        if options['filepath'] == None:
            raise FileNotFoundError

        map = pickle.load(open(options['filepath'], 'rb'))

        return map