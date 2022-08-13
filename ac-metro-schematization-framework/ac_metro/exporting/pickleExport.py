from  ac_metro.exporting.abstractExport import AbstractExport
import pickle

class PickleExport(AbstractExport):
    '''
    Load a metro map from an Graphml file.
    '''

    def export(map, options:dict):
        if options['filepath'] == None:
            raise FileNotFoundError

        pickle.dump(map, open(options['filepath'], 'wb'))
