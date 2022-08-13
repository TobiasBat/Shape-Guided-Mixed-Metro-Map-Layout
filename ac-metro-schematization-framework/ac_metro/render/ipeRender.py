from ac_metro.error.metroErrors import OptionsError
from ac_metro.render.abstractRender import AbstractRender
from ac_metro.map.map import Map
from ipey.document import Document
from ipey.primitive import Line, Glyph
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


def render_map(map, document):
    page = document.createPage()
    
    for key, stations in map.getLines():
        print(f"line key is {key}")
        print(f"map line colors are {map.lineColors}")
        points = stations_to_coords(map, stations)
        if points == -1:
            continue
        line = Line(stations_to_coords(map, stations))
        line.stroke = (map.lineColors[key] if key in map.lineColors else '#000000')
        line.pen = "4"
        line.layer = f'{key}'
        page.add(line)

    for v, data in map.graph.nodes(data=True):
        station = data['station']
        glyph = Glyph((station.x, station.y), type='mark/fpill_1(sfx)')
        glyph.stroke = "#000000"
        glyph.fill = "#ffffff"
        glyph.size = "6" #"0.6"
        glyph.layer = 'Station'
        page.add(glyph)

class IpeRender(AbstractRender):
    '''
    Use IPE to render the map.
    '''

    def render(map: Map, options=None):
        '''
        
        '''

        if not options:
            raise OptionsError("IpeRender: No options given")

        if not options['filepath']:
            raise OptionsError("IpeRender: No path in options given")

        if not options['filename']:
            raise OptionsError("IpeRender: No filename in options given")

        document = Document(styles=["ac_metro/render/static/metro_marks.isy"])
        document.crop = True

        print(f"called ipe render with name {options['filename']}")

        if 'render_all_snapshots' in options and options['render_all_snapshots'] == True:
            for snapshot in map.snapshots:
                ColorDictUtility.execute(snapshot.map, {"preset" : options["filename"]})
                render_map(snapshot.map, document)
        else:
            render_map(map, document)

        ipename = f'{options["filepath"]}/{options["filename"]}.ipe'
        print(f"writing to {ipename}")
        document.write(ipename)
        return 0