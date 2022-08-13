from ac_metro.utility.abstractUtility import AbstractUtility

class ColorDictUtility(AbstractUtility):
    '''
    Use a dictionary with identifier-color pairs to color a metro map.
    '''
    def execute(map, options=None):
        if options == None:
            return

        print(options)

        if "preset" in options:
            if options["preset"] == "Tiny" or options["preset"] == "Micro":
                map.lineColors = {"l0" : "#46b446", "l1" : "#fabe14"}
            elif options["preset"] == "Vienna":
                map.lineColors = {"U1" : "#D52230", "U2" : "#935F96", "U3" : "#E47835", "U4" : "#16984C", "U6" : "#8A633F"}
            elif options["preset"] == "Paris" or options["preset"] == "ParisE" or options["preset"] == "ParisCl" or options["preset"] == "ParisCi":
                map.lineColors = {"l0": "#f4d33a", "l1": "#2625d8", "l2": "#ff8000", "l3": "#87d3df", "l4": "#c24093",
                                  "l5": "#f07f4b", "l6": "#78bc8b", "l7": "#cab754", "l8": "#ff6a46", "l9": "#ceb325",
                                  "l10": "#06a02c", "l11": "#6dc5e0", "l12": "#6dc5e0", "l13": "#6c1786",
                                  "l14": "#ef9fb5", "l15": "#ef9fb5", "l16": "#78cd73", "l17": "#8e5b3c",
                                  "l18": "#c2a1ca", "l19": "#000000", "l20": "#000000"}
            elif options["preset"] == "Montreal":
                map.lineColors = {"l0" : "#D52230", "l1" : "#935F96", "l2" : "#E47835", "l3" : "#16984C"}
            elif options["preset"] == "Prague":
                map.lineColors = {"l0" : "#46b446", "l1" : "#fabe14", "l2" : "#f51923"}
            elif options["preset"] == "Taipei":
                map.lineColors = {"l0" : "#a74c00", "l1" : "#ff0000", "l2" : "#f48b9a", "l3" : "#008000", "l4" : "#b5d931", "l5" : "#ffa500", "l6" : "#0018ff", "l7" : "#ffa500", "l8" : "#ffa500", "l9" : "#d1e234"}
            elif options["preset"] == "Moscow":
                map.lineColors = {"l0": "#e82f2f", "l1": "#3db556", "l2": "#0076c0", "l3": "#29c4f0", "l4": "#804000", "l5": "#ff801f", "l6": "#853892", "l7": "#efc000", "l8": "#e5b600", "l9": "#9f9f9a", "l10": "#9fc000", "l11": "#6fc0c7", "l12": "#88acdd", "l13": "#477cb2", "l14": "#cc9966", "l15": "#e384ce"}
            elif options["preset"] == "Lissabon" or options["preset"] == "Lisbon":
                map.lineColors = {"l0" : "#D4396C", "l1" : "#F0BC41", "l2" : "#43989C", "l3" : "#5D84BF", "l4" : "#C8C8B4", "l5" : "#C8C8B4", "l6" : "#C8C8B4", "l7" : "#C8C8B4", "l8" : "#C8C8B4"}
            elif options["preset"] == "Berlin" or options["preset"] == "BerlinSU" or options["preset"] == "BerlinS":
                map.lineColors = {"l0": "#54a721", "l1": "#ff3700", "l2": "#019276", "l3": "#ffd900", "l4": "#662e16", "l5": "#6e4d9b", "l6": "#358fbf", "l7": "#0a3b84", "l8": "#ff7200", "l9": "#e6a0ff", "l10": "#00501e", "l11": "#146432", "l12": "#146432", "l13": "#4b50e1", "l14": "#913705", "l15": "#c8a06e", "l16": "#c8a06e", "l17": "#c8a06e", "l18": "#ff8400", "l19": "#824be1", "l20": "#824be1", "l21": "#4be164", "l22": "#820505", "l23": "#000000", "l24": "#000000"}
            # elif options["preset"] == "BerlinS":
            #     map.lineColors = {"l0" : "#9696FF"}
            elif options["preset"] == "Tokyo" or options["preset"] == "TokyoCut":
                map.lineColors = {"l0": "#ff0000", "l1": "#ff0000", "l2": "#ff9900", "l3": "#cccccc", "l4": "#6699ff",
                                  "l5": "#006600", "l6": "#000000", "l7": "#000000", "l8": "#000000", "l9": "#ffcc00",
                                  "l10": "#000000", "l11": "#000000", "l12": "#9999ff", "l13": "#000000",
                                  "l14": "#339999", "l15": "#996633", "l16": "#ee5588", "l17": "#0077bb",
                                  "l18": "#66bb55", "l19": "#cc0077", "l20": "#cc0077", "l21": "#000000",
                                  "l22": "#000000", "l23": "#000000"}
            elif options["preset"] == "Sydney":
                map.lineColors = { "l1" : "#ffd700", "l2" : "#ff0000", "l3" : "#000080", "l4" : "#ff8c00", "l5" : "#228b22", "l6" : "#6495ed", "l7" : "#6a5acd", "l8" : "#c71585", "l9" : "#87cefa", "l10" : "#000000"}
            elif options["preset"] == "Karlsruhe":
                map.lineColors = {"l0": "#ee1c23", "l1": "#0072bc", "l2": "#937138", "l3": "#fecb09", "l4": "#16c1f3", "l5": "#7fc241", "l11": "#f7931d", "l6": "#00a76e", "l7": "#a066ab", "l9": "#9e174c", "l10": "#f69795"}
            elif options["preset"] == "KA":
                map.lineColors = {"KA_S1": "#203f6b","KA_S2": "#96a44e","KA_S3": "#da913e","KA_S4": "#b82728","KA_S5": "#edd64a"}
            elif options["preset"] == "VA":
                map.lineColors = {"VA_S1": "#992426","VA_S2": "#614967","VA_S3": "#a8b861","VA_S4": "#f0bc19","VA_S5": "#e6b5c8","VA_S7": "#3b7fbc"}
            elif options["preset"] == "TI":
                map.lineColors = {"TI_S2": "#dca9bd","TI_S3": "#afafaf","TI_S4": "#65476d","TI_S5": "#64b4a6","TI_S6": "#9d3541","TI_S7": "#ab919c","TI_S8": "#40add6"}
            elif options["preset"] == "SA":
                map.lineColors = {"SA_S1": "#c13c4b","SA_S11": "#c13c4b","SA_S2": "#246999","SA_S3": "#94a837","SA_S4": "#6f4479",}
            elif options["preset"] == "OO":
                map.lineColors = {"OO_S1": "#da9c3b","OO_S2": "#5fa092","OO_S3": "#405b7a","OO_S4": "#a1b63d","OO_S5": "#eab1c9"}
            elif options["preset"] == "ST":
                map.lineColors = {"ST_S1": "#a3b34d","ST_S11": "#a3b34d","ST_S3": "#ddadc2","ST_S31": "#ddadc2","ST_S5/S51": "#5e4268","ST_S5": "#5e4268","ST_S51": "#5e4268","ST_S6": "#e09231","ST_S61": "#e09231","ST_S7": "#942936","ST_S8": "#65b8a4","ST_S9": "#b9aeba"}
            elif options["preset"] == "WI":
                map.lineColors = {"WI_S1": "#dfa7bf","WI_S2": "#dfa7bf","WI_S3": "#dfa7bf","WI_S4": "#dfa7bf","WI_S7": "#dfa7bf","WI_S40": "#4ab7e9","WI_S45": "#9db53a","WI_S50": "#4ab7e9","WI_S60": "#4ab7e9","WI_S80": "#4ab7e9"}
            elif options["preset"] == "AT_S_Bahn":
                map.lineColors = {"KA_S1": "#203f6b","KA_S2": "#96a44e","KA_S3": "#da913e","KA_S4": "#b82728","KA_S5": "#edd64a","VA_S1": "#992426","VA_S2": "#614967","VA_S3": "#a8b861","VA_S4": "#f0bc19","VA_S5": "#e6b5c8","VA_S7": "#3b7fbc","TI_S2": "#dca9bd","TI_S3": "#afafaf","TI_S4": "#65476d","TI_S5": "#64b4a6","TI_S6": "#9d3541","TI_S7": "#ab919c","TI_S8": "#40add6","SA_S1": "#c13c4b","SA_S11": "#c13c4b","SA_S2": "#246999","SA_S3": "#94a837","SA_S4": "#6f4479","OO_S1": "#da9c3b","OO_S2": "#5fa092","OO_S3": "#405b7a","OO_S4": "#a1b63d","OO_S5": "#eab1c9","ST_S1": "#a3b34d","ST_S11": "#a3b34d","ST_S3": "#ddadc2","ST_S31": "#ddadc2","ST_S5": "#5e4268","ST_S51": "#5e4268","ST_S6": "#e09231","ST_S61": "#e09231","ST_S7": "#942936","ST_S8": "#65b8a4","ST_S9": "#b9aeba","WI_S1": "#dfa7bf","WI_S2": "#dfa7bf","WI_S3": "#dfa7bf","WI_S4": "#dfa7bf","WI_S7": "#dfa7bf","WI_S40": "#4ab7e9","WI_S45": "#9db53a","WI_S50": "#4ab7e9","WI_S60": "#4ab7e9","WI_S80": "#4ab7e9"}
            elif options["preset"] == "Singapore":
                map.lineColors = {"l0": "#1a9f39", "l1": "#d82e29", "l2": "#1a9f39", "l3": "#8c5091", "l4": "#f4b135", "l5": "#f4b135", "l6": "#824b0a", "l7": "#0000ff", "l8": "#0000ff"}
        else:
            for key, value in options.items():
                map.lineColors[key] = value
