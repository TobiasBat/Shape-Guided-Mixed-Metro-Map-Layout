# Helper methods and constants
import math

import numpy as np

from ac_metro.importing.graphmlImport import GraphmlImport
from ac_metro.render.matplotlibGridRender import MatplotlibGridRender

I = u'\u2502'
L = u'\u2514'
D = u'\u2500'
E = u'\u257C'
LINE = L + D + D + D + E + ' '
HOOK = I + '\n' + LINE
S_HOOK = '      ' + I + '\n      ' + LINE
S_HOOK_S = '      ' + LINE

def set_meta_data3(sys_args):
    name = None
    directory = None

    if sys_args[1] == "B":
        name = "BerlinS"
        directory = "berlinS"
    elif sys_args[1] == "BU":
        name = "BerlinSU"
        directory = "berlinSU"
    elif sys_args[1] == "BUC":
        name = "BerlinSU-Cirlce"
        directory = "berlinSU-Circle"
    elif sys_args[1] == "BB":
        name = "Berlin"
        directory = "berlinS-Bear"
    elif sys_args[1] == "BUB":
        name = "Berlin"
        directory = "berlinSU-Bear"
    elif sys_args[1] == "T":
        name = "Taipei"
        directory = "taipei"
    elif sys_args[1] == "To":
        name = "Tokyo"
        directory = "tokyo"
    elif sys_args[1] == "L":
        name = "Lisbon"
        directory = "lisbon"
    elif sys_args[1] == "M":
        name = "Moscow"
        directory = "moscow"
    elif sys_args[1] == "Mo":
        name = "Montreal"
        directory = "montreal"
    elif sys_args[1] == "Pr":
        name = "Prague"
        directory = "prague"
    elif sys_args[1] == "Pa":
        name = "Paris"
        directory = "paris"
    elif sys_args[1] == "PaCi":
        name = "ParisCi"
        directory = "paris-Circle"
    elif sys_args[1] == "PaCl":
        name = "ParisCl"
        directory = "paris-Cloud"
    elif sys_args[1] == "PaE":
        name = "ParisE"
        directory = "paris-Eye"
    elif sys_args[1] == "S-C":
        name = "Singapore"
        directory = "singapore-Circle"
    elif sys_args[1] == "S-H":
        name = "Singapore"
        directory = "singapore-Heart"

    return name, directory

def set_meta_data2(sys_args):
    name = None
    directory = None

    if sys_args[1] == "B":
        name = "BerlinS"
        directory = "BerlinS"
    elif sys_args[1] == "Bg":
        name = "BerlinS"
        directory = "BerlinS-geo"
    elif sys_args[1] == "Bs":
        name = "BerlinS"
        directory = "BerlinS-smooth"
    elif sys_args[1] == "BU":
        name = "BerlinSU"
        directory = "BerlinSU"
    elif sys_args[1] == "BUg":
        name = "BerlinSU"
        directory = "BerlinSU-geo"
    elif sys_args[1] == "BUs":
        name = "BerlinSU"
        directory = "BerlinSU-smooth"
    elif sys_args[1] == "BUC":
        name = "BerlinSU-Cirlce"
        directory = "BerlinSU-Circle"
    elif sys_args[1] == "BB":
        name = "Berlin"
        directory = "BerlinS-Bear"
    elif sys_args[1] == "BUB":
        name = "Berlin"
        directory = "BerlinSU-Bear"
    elif sys_args[1] == "T":
        name = "Taipei"
        directory = "Taipei"
    elif sys_args[1] == "Tg":
        name = "Taipei"
        directory = "Taipei-geo"
    elif sys_args[1] == "Ts":
        name = "Taipei"
        directory = "Taipei-smooth"
    elif sys_args[1] == "Tns":
        name = "Taipei"
        directory = "Taipei-no-shape"
    elif sys_args[1] == "To":
        name = "Tokyo"
        directory = "Tokyo"
    elif sys_args[1] == "Tos":
        name = "Tokyo"
        directory = "Tokyo-smooth"
    elif sys_args[1] == "Tog":
        name = "Tokyo"
        directory = "Tokyo-geo"
    elif sys_args[1] == "ToC":
        name = "Tokyo"
        directory = "TokyoCut"
    elif sys_args[1] == "ToS":
        name = "TokyoSmall"
        directory = "TokyoSmall"
    elif sys_args[1] == "L":
        name = "Lisbon"
        directory = "Lisbon"
    elif sys_args[1] == "Lg":
        name = "Lisbon"
        directory = "Lisbon-geo"
    elif sys_args[1] == "Ls":
        name = "Lisbon"
        directory = "Lisbon-smooth"
    elif sys_args[1] == "M":
        name = "Moscow"
        directory = "Moscow"
    elif sys_args[1] == "Mg":
        name = "Moscow"
        directory = "Moscow-geo"
    elif sys_args[1] == "Ms":
        name = "Moscow"
        directory = "Moscow-smooth"
    elif sys_args[1] == "Mo":
        name = "Montreal"
        directory = "Montreal"
    elif sys_args[1] == "Mog":
        name = "Montreal"
        directory = "Montreal-geo"
    elif sys_args[1] == "Mos":
        name = "Montreal"
        directory = "Montreal-smooth"
    elif sys_args[1] == "Mo":
        name = "Montreal"
        directory = "Montreal"
    elif sys_args[1] == "Pr":
        name = "Prague"
        directory = "Prague"
    elif sys_args[1] == "Pa":
        name = "Paris"
        directory = "Paris"
    elif sys_args[1] == "Pag":
        name = "Paris"
        directory = "Paris-geo"
    elif sys_args[1] == "PaCi":
        name = "ParisCi"
        directory = "Paris-Circle"
    elif sys_args[1] == "PaCis":
        name = "ParisCi"
        directory = "Paris-Circle-smooth"
    elif sys_args[1] == "PaCl":
        name = "ParisCl"
        directory = "Paris-Cloud"
    elif sys_args[1] == "PaE":
        name = "ParisE"
        directory = "Paris-Eye"
    elif sys_args[1] == "PaEs":
        name = "ParisE"
        directory = "Paris-Eye-smooth"
    elif sys_args[1] == "Sinus":
        name = "Sinus"
        directory = "Sin"
    elif sys_args[1] == "S-C":
        name = "Singapore"
        directory = "Singapore-Circle"
    elif sys_args[1] == "S-H":
        name = "Singapore"
        directory = "Singapore-Heart"
    elif sys_args[1] == "S-Hs":
        name = "Singapore"
        directory = "Singapore-Heart-smooth"
    elif sys_args[1] == "Sg":
        name = "Singapore"
        directory = "Singapore-geo"
    elif sys_args[1] == "Ss":
        name = "Singapore"
        directory = "Singapore-smooth"

    return name, directory

def set_meta_data(sys_args):
    name = None
    directory = None
    ins_tech = None
    input_type = None

    if sys_args[1] == "B":
        name = "Berlin"
        directory = "Berlin"
        ins_tech = "distance"
    elif sys_args[1] == "BPB":
        name = "Berlin"
        directory = "Berlin-Planar-Bear"
        ins_tech = "distance"
    elif sys_args[1] == "BS":
        name = "BerlinS"
        directory = "BerlinS"
        ins_tech = "distance"
    elif sys_args[1] == "BSUNP":
        name = "Berlin"
        directory = "BerlinS+U-NotPlanar"
        ins_tech = "distance"
    elif sys_args[1] == "BSUP":
        name = "Berlin"
        directory = "BerlinS+U-Planar"
        ins_tech = "distance"
    elif sys_args[1] == "L":
        name = "Lissabon"
        directory = "Lissabon"
        ins_tech = "distance"
    elif sys_args[1] == "M":
        name = "Moscow"
        directory = "Moscow-Planar"
        ins_tech = "distance"
    elif sys_args[1] == "T":
        name = "Taipei"
        directory = "Taipei"
        ins_tech = "distance"
    elif sys_args[1] == "To":
        name = "Tokyo"
        directory = "Tokyo"
        ins_tech = "distance"
    elif sys_args[1] == "ToC":
        name = "Tokyo"
        directory = "TokyoCut"
        ins_tech = "distance"
    elif sys_args[1] == "ToS":
        name = "TokyoSmall"
        directory = "TokyoSmall"
        ins_tech = "distance"
    elif sys_args[1] == "Tiny":
        name = "Tiny"
        directory = "Tiny"
        # ins_tech = "crossings"
        ins_tech = "distance"
    elif sys_args[1] == "Micro":
        name = "Micro"
    else:
        print("No Input defined")
        quit()
    if sys_args[2] == "g":
        input_type = "geo"
    elif sys_args[2] == "s":
        input_type = "smooth"
    elif sys_args[2] == "m":
        input_type = "mixed"
    else:
        print("No Input Type defined")
        quit()
    return name, directory, ins_tech, input_type


def is_smooth(edge):
    return edge.hasContextualInformation(GraphmlImport.IMPORT_NAME()) \
           and "smooth" in edge.getContextualInformation(GraphmlImport.IMPORT_NAME()) \
           and edge.getContextualInformation(GraphmlImport.IMPORT_NAME())["smooth"]

def is_turning(station):
    return station.hasContextualInformation(GraphmlImport.IMPORT_NAME()) \
           and "turningpoint" in station.getContextualInformation(GraphmlImport.IMPORT_NAME()) \
           and station.getContextualInformation(GraphmlImport.IMPORT_NAME())["turningpoint"]

def is_smooth_station(station):
    return station.hasContextualInformation(GraphmlImport.IMPORT_NAME()) \
           and "smoothNode" in station.getContextualInformation(GraphmlImport.IMPORT_NAME()) \
           and station.getContextualInformation(GraphmlImport.IMPORT_NAME())["smoothNode"]

def dist(x1, y1, x2, y2):
    '''
    Returns the Euclidean distance
    '''
    return math.sqrt(math.pow(x1 - x2, 2) + math.pow(y1 - y2, 2))


def compute_x_axis_angle(x1, y1, x2, y2):
    '''
    Returns the angle between the x axis and a line throught the points
    '''
    return math.atan2(y2 - y1, x2 - x1)


def compute_bend_angle(x1, y1, x2, y2, x3, y3):
    '''
    computes angle between lines through point 1 and 2 and through points 2 and 3
    '''
    predCoord = np.array([x1, y1])
    currCoord = np.array([x2, y2])
    nextCoord = np.array([x3, y3])

    # print(f"Angle: {(round(x1, 1), round(y1, 1))} - {(round(x2, 1), round(y2, 1))} - {(round(x3, 1), round(y3, 1))}")

    uv = predCoord - currCoord
    vw = nextCoord - currCoord
    uvD = np.linalg.norm(uv)
    vwD = np.linalg.norm(vw)
    uv = uv / uvD
    vw = vw / vwD

    dot = round(np.dot(uv, vw), ndigits=6)

    # print(f"is {np.arccos(dot)}\n")
    # input()
    return np.arccos(dot)


def angle_cost(angle):
    '''
    Different function, which assign a penalty value to an angle. Used for Port version.
    '''
    # TODO CHECK MAGIC CONSTANT!
    return math.pow(math.pi - angle, 2) + 20


def angle_cost_2(angle):
    '''
    Different function, which assign a penalty value to an angle. Used for Port version.
    '''
    # math.pi is the maximal angle between two lines (and thereby 0 is the minimal angle, which gets a penalty of math.pi)
    # so this is simply 2*math.pi - angle
    return (2*math.pi - angle) #+ math.pi


def bprint(text, end="\n"):
    print("\033[1m" + str(text) + "\033[0m", end=end)


def render_grid(grid, path, dpi=300):
    MatplotlibGridRender.render(None, {
        "grid": grid,
        "just_grid": True,
        "grid_labels": True,
        "show_nodes": False,
        "font_size": 2,
        "dpi": dpi,
        # "grid_edge_labels": "weight",
        "to_file": path
    })


def render_map(map, grid, path):
    MatplotlibGridRender.render(map, {
        "grid": grid,
        "just_grid": False,
        "grid_labels": False,
        "show_nodes": False,
        "station_names": True,
        "to_file": path
    })


# Mark nodes, which are too close as not valid for candidate sets. DOESNT SEEM TO WORK RIGHT
def mark_not_valid_candidates(grid, treshhold):
    marked = set()
    for node, data in grid.nodes(data=True):
        if not data["shape_part"]:
            continue
        if node in marked:
            continue
        x1, y1 = data["x"], data["y"]
        for node2, data2 in grid.nodes(data=True):
            if node2 in marked or node == node2:
                continue
            x2, y2 = data2["x"], data2["y"]
            if dist(x1, y1, x2, y2) < treshhold:
                grid.nodes()[node2]["not_candidate"] = True
                marked.add(node2)
