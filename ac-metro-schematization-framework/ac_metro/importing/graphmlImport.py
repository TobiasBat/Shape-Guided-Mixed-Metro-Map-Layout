import csv
import networkx as nx
from  ac_metro.importing.abstractImport import AbstractImport
from ac_metro.map.map import Map, Station

class GraphmlImport(AbstractImport):
    '''
    Load a metro map from an Graphml file.
    '''
    noLineIDCounter = 1

    @staticmethod
    def IMPORT_NAME():
        return 'GraphmlImport'

    def load(options:dict):
        if options['filepath'] == None:
            raise FileNotFoundError

        G_in = nx.read_graphml(options['filepath'])

        if 'x_attr' in options:
            xAttr = options['x_attr']
        else:
            xAttr = 'x'

        if 'y_attr' in options:
            yAttr = options['y_attr']
        else:
            yAttr = 'y'

        idAttr = None
        if 'id_attr' in options:
            idAttr = options['id_attr']

        labelAttr = None
        if 'label_attr' in options:
            labelAttr = options['label_attr']
      
        if 'line_attr' in options:
            lineAttr = options['line_attr']
        else:
            lineAttr = 'line'

        contextual_station_information_map = {}
        if 'contextual_station_information' in options:
            contextual_station_information_map = options['contextual_station_information']

        contextual_edge_information_map = {}
        if 'contextual_edge_information' in options:
            contextual_edge_information_map = options['contextual_edge_information']

        map = Map()

        for v, data in G_in.nodes(data=True):
            x = float(data[xAttr])
            y = float(data[yAttr])
            
            if idAttr:
                id = data[idAttr]
            else:
                id = v

            if labelAttr:
                lbl = data[labelAttr]

            contextual_station_information = {}
            for contextual_station_information_name, contextual_station_information_attr in contextual_station_information_map.items():
                contextual_station_information[contextual_station_information_name] = data[contextual_station_information_attr]

            s = Station(id, x=x, y=y, label=lbl, author=GraphmlImport.IMPORT_NAME(), contextualInformation=contextual_station_information)

            map.addStation(s)

        for u, v, data in G_in.edges(data=True):
            contextual_edge_information = {}
            for contextual_edge_information_name, contextual_edge_information_attr in contextual_edge_information_map.items():
                contextual_edge_information[contextual_edge_information_name] = data[contextual_edge_information_attr]

            if lineAttr in data:
                for line in data[lineAttr].split(','):
                    # print(f"Adding {u} - {v} to line {line.strip()}")
                    map.addConnection(u, v, line.strip(), GraphmlImport.IMPORT_NAME(), contextual_edge_information)
            else:
                print(f"could not find a line for {u} - {v}")
                map.addConnection(u,v, f'NL_{GraphmlImport.noLineIDCounter}', GraphmlImport.IMPORT_NAME(), contextual_edge_information)
                GraphmlImport.noLineIDCounter += 1

        return map
