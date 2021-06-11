
def add_3d_lines(self, view, i, j, color, viewer=None):
    xi, yi, zi = self.get_xyz_coord(i)
    xj, yj, zj = self.get_xyz_coord(j)
    cylinder_specs = {"start": {"x": xi, "y": yi, "z": zi},
                      "end":  {"x": xj, "y": yj, "z": zj},
                      "radius": 0.5,
                      "fromCap": 2,
                      "toCap": 2,
                      "color": color}
    if viewer is not None:
        view.addCylinder(cylinder_specs, viewer=viewer)
    else:
        view.addCylinder(cylinder_specs)


def plot_3d_data(self, view, ij_data, metric=None, viewer=None,
                 **kwargs):
    self.filter_ij_data(ij_data, "pdb", **kwargs)
    ij_colors = self.get_ij_colors(ij_data, metric)
    for i, j, color in zip(*ij_colors):
        color = mp.colors.rgb2hex(color)
        color = "0x"+color[1:]
        window = self.window[ij_data]
        for w in range(window):
            io = i+w
            jo = j+window-1-w
            if io in self.pdb_validres and jo in self.pdb_validres:
                self.add_3d_lines(view, io, jo, color, viewer)


def set_3d_colors(self, view, colorby, viewer=None):
    colorby_function = getattr(self, "get_colorby_"+colorby)
    colors = colorby_function("pdb")
    color_selector = {}
    for res in self.pdb.get_residues():
        res = res.get_id()
        color = colors[res[1]-1]
        if color in color_selector.keys():
            color_selector[color].append(res[1])
        else:
            color_selector[color] = [res[1]]
    for color in color_selector.keys():
        selector = {'resi': color_selector[color]}
        style = {"cartoon": {"color": color, "opacity": 0.8}}
        if viewer is None:
            view.setStyle(selector, style)
        else:
            view.setStyle(selector, style, viewer=viewer)
    selector = {'resi': self.pdb_validres, 'invert': 'true'}
    style = {"cross": {"hidden": "true"}}
    if viewer is None:
        view.setStyle(selector, style)
    else:
        view.setStyle(selector, style, viewer=viewer)


def set_3d_view(self, view):
    with open(self.paths["pdb"], 'r') as pdb_file:
        pdb_str = pdb_file.read()
    view.addModel(pdb_str, 'pdb')
    view.setStyle({'chain': 'A'}, {"cartoon": {'color': 'grey'}})
    view.zoomTo()
    return view


def make_3d(self, view=None, viewer=None, ij_data=None, metric=None,
            colorby="sequence", **kwargs):
    if view is None:
        view = py3Dmol.view()
        view = self.set_3d_view(view)
    if hasattr(self, "profile"):
        colorby = "profile"
    self.set_3d_colors(view, colorby, viewer)
    if metric == "Distance":
        self.set_3d_distances(ij_data)
    if ij_data is not None:
        self.plot_3d_data(view, ij_data, metric, viewer, **kwargs)
    return view
