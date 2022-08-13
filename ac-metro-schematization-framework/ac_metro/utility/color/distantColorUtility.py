from ac_metro.utility.abstractUtility import AbstractUtility

from colormath.color_objects import LabColor, sRGBColor
from colormath.color_diff import delta_e_cie1976
from colormath.color_conversions import convert_color
import networkx as nx

COLORS = ['#4E79A7', '#E15759', '#59A14F', '#D37295', '#9D7660', '#F28E2B', '#499894', '#79706E', '#B07AA1', '#B6992D', '#A0CBE8', '#FFBE7D', '#8CD17D', '#F1CE63', '#86BCB6', '#FF9D9A', '#BAB0AC', '#FABFD2', '#D4A6C8', '#D7B5A6']
NUM_C = len(COLORS)

class DistantColorUtility(AbstractUtility):
    '''
    Color a map by a greedy approach that assigns lines a color such that the perceptional difference is maximized. 
    '''

    def execute(map, options=None):
        cGraph = nx.Graph()

        COLORS_LAB = []
        C_DIFF = {}
        for color in COLORS:
            rgb = sRGBColor.new_from_rgb_hex(color)
            lab = convert_color(rgb, LabColor)
            COLORS_LAB.append(lab)

        for i in range(0, len(COLORS)):
            for j in range(i + 1, len(COLORS)):
                c1 = COLORS[i]
                c2 = COLORS[j]
                l1 = COLORS_LAB[i]
                l2 = COLORS_LAB[j]

                diff = delta_e_cie1976(l1,l2)

                C_DIFF[(c1,c2)] = diff
                C_DIFF[(c2,c1)] = diff
                COLORS_LAB.append(lab)

        num = len(map.lines)
        repeat = int(num / NUM_C) + 1

        colors = []
        colorsLab = []
        for i in range(repeat):
            colors.extend(COLORS)
            colorsLab.extend(COLORS_LAB)

        colors = colors[0:num]
        colorsLab = colorsLab[0:num]
        assignedColors = []

        intersections = {}
        for key1, line1 in map.getLines():
            for key2, line2 in map.getLines():

                if key1 == key2:
                    continue

                valIntersection = DistantColorUtility.intersection(line1, line2)

                intersections[(key1,key2)] = valIntersection
                intersections[(key2,key1)] = valIntersection

                cGraph.add_edge(key1, key2, intersection = valIntersection)

        #print(cGraph.nodes())
        #print(cGraph.edges())
        
        edges = sorted(cGraph.edges(data=True), key=lambda t: t[2].get('intersection', 1))

        for (key1, key2, value) in reversed(edges):
            #print(key1, key2, value)
            hasColor1 = False
            hasColor2 = False

            if 'color' in cGraph.nodes[key1]:
                hasColor1 = True

            if 'color' in cGraph.nodes[key2]:
                hasColor2 = True

            if hasColor1 and hasColor2:
                continue

            total1 = cGraph.degree(key1, weight='intersection')
            total2 = cGraph.degree(key2, weight='intersection')

            #process higher degree first
            if total2 > total1:
                temp = key1
                key1 = key2
                key2 = temp
                temp = hasColor1
                hasColor1 = hasColor2
                hasColor2 = temp

            if not hasColor1:

                colDiff = [0] * len(colors)
                for n in cGraph.neighbors(key1):
                    weight = cGraph.edges[(key1,key2)]['intersection']
                    #print(key1, n)
                    if 'color' in cGraph.nodes[n]:
                        colorN = cGraph.nodes[n]['colorLab']
                        #print('hasColor')
                        for index,col in enumerate(colorsLab):
                            colDiff[index] += delta_e_cie1976(colorN, col) * weight

                if max(colDiff) <= 0:
                    for index,col1 in enumerate(colorsLab):
                        for col2 in assignedColors:
                            colDiff[index] -= delta_e_cie1976(col1, col2)

                
                ind = colDiff.index(max(colDiff))
                cGraph.nodes[key1]['color'] = colors[ind]
                cGraph.nodes[key1]['colorLab'] = colorsLab[ind]

                assignedColors.append(colorsLab[ind])
                #print(colors[ind])

                colors.pop(ind)
                colorsLab.pop(ind)

            if not hasColor2:

                colDiff = [0] * len(colors)
                for n in cGraph.neighbors(key2):
                    weight = cGraph.edges[(key1,key2)]['intersection']

                    if 'color' in cGraph.nodes[n]:
                        colorN = cGraph.nodes[n]['colorLab']

                        for index,col in enumerate(colorsLab):
                            colDiff[index] += delta_e_cie1976(colorN, col) * weight

                if max(colDiff) <= 0:
                    for index,col1 in enumerate(colorsLab):
                        for col2 in assignedColors:
                            colDiff[index] -= delta_e_cie1976(col1, col2)

                #print(max(colDiff))

                ind = colDiff.index(max(colDiff))
                cGraph.nodes[key2]['color'] = colors[ind]
                cGraph.nodes[key2]['colorLab'] = colorsLab[ind]

                assignedColors.append(colorsLab[ind])
                #print(colors[ind])
                colors.pop(ind)
                colorsLab.pop(ind)

        colorDict = {}
        for key, line in map.getLines():
            map.lineColors[key] = cGraph.nodes[key]['color']


    @staticmethod
    def intersection(lst1, lst2):
        return len(list(set(lst1) & set(lst2))) 