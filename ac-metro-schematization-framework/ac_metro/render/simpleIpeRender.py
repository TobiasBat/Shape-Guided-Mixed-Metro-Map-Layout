from ac_metro.error.metroErrors import OptionsError
from ac_metro.render.abstractRender import AbstractRender
from ac_metro.map.map import Map
from ipey.document import Document
from ipey.primitive import Line, Glyph, Label
from ac_metro.utility.color.colorDictUtility import ColorDictUtility


def stations_to_coords(map, stations):
    if len(stations) < 2:
        return -1
    last = stations[0]
    st_last = map.graph.nodes[last]['station']
    coords = [(st_last.x, st_last.y)]
    for s in stations[1:]:
        st_new = map.graph.nodes[s]['station']
        control_points = map.graph[last][s]["edge"].getControlPoints(last)
        if len(control_points) == 0:
            coords.append((st_new.x, st_new.y))
        else:
            coords.append((control_points[0][0], control_points[0][1]))
            for point in control_points[1:]:
                coords.append((point[0], point[1]))
            coords.append((st_new.x, st_new.y))
        last = s

    if map.graph.has_edge(last, stations[0]):
        st_new = map.graph.nodes[stations[0]]['station']
        control_points = map.graph[last][stations[0]]["edge"].getControlPoints(last)
        if len(control_points) == 0:
            coords.append((st_new.x, st_new.y))
        else:
            coords.append((control_points[0][0], control_points[0][1]))
            for point in control_points[1:]:
                coords.append((point[0], point[1]))
                coords.append((st_new.x, st_new.y))
    return coords


def render_map(map, document, options):
    page = document.createPage()

    print(f"map has the lines:\n{map.getLines()}")
    
    for key, stations in map.getLines():
        # print(f"line key is {key} with {stations}")
        if key.startswith("deg2heu_internal_single_stop_helper_line_"):
            # print('skipping')
            continue
        # print(f"coloring {key} which is in map.lineColors: {key in map.lineColors}, i.e., {map.lineColors[key]}")
        line = Line(stations_to_coords(map, stations))
        line.stroke = map.lineColors[key]
        # print(map.lineColors)
        # print(key)
        line.pen = "1"
        line.layer = f'{key}'
        page.add(line)
        print(f"added a line {key}")

    print(page.Scene)
    page.createLayer("Labels")
    page.createLayer("Glyphs")
    for v, data in map.graph.nodes(data=True):
        station = data['station']
        if "use_labels" in options and options["use_labels"]:
            text = station.label
            text = text.replace("_", "-")
            text = text.replace("\\", "")
            text = text.replace("Ã", "A")
            text = text.replace(u"\u00c3", "-00c3-")
            text = text.replace(u"\u2022", "-2022-")
            text = text.replace(u"\u0083", "-0083-")
            text = text.replace(u"\u008A", "-008A-")
            text = text.replace(u"\u0082", "-0082-")
            text = text.replace(u"\u0087", "-0087-")
            text = text.replace(u"\u0093", "-0093-")
            text = text.replace(u"\u0089", "-0089-")
            label = Label(text, (station.x, station.y))
            label.layer = "Labels"
            page.add(label)
        glyph = Glyph((station.x, station.y))
        glyph.layer = "Glyphs"
        page.add(glyph)

class SimpleIpeRender(AbstractRender):
    '''
    Use IPE to render the map.
    '''

    @staticmethod
    def render(metro_map: Map, options=None):
        '''
        
        '''

        if not options:
            raise OptionsError("IpeRender: No options given")

        if not options['filepath']:
            raise OptionsError("IpeRender: No path in options given")

        if not options['filename']:
            raise OptionsError("IpeRender: No filename in options given")

        document = Document()
        document.crop = True

        print(f"called simple ipe render with name {options['filename']}")

        ColorDictUtility.execute(metro_map, {"preset": options["preset"] if "preset" in options else options["filename"]})

        render_map(metro_map, document, options)

        ipename = f'{options["filepath"]}/{options["filename"]}.ipe'
        print(f"writing to {ipename}")
        document.write(ipename)