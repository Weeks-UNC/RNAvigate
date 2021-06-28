from plots import get_rows_columns, view_colormap
import py3Dmol
import matplotlib.colors as mpc


class Mol():
    def __init__(self, pdb, profiles=[], ijs=[], labels=[]):
        self.pdb = pdb
        self.labels = []
        self.profiles = []
        self.ij_windows = []
        self.ijs = []
        for sample in zip(ijs, profiles, labels):
            self.add_sample(sample)

    def add_sample(self, ij, profile=None, label=None):
        self.ijs.append(ij.get_ij_colors())
        self.ij_windows.append(ij.window)
        view_colormap(ij)
        if label is not None:
            self.labels.append(label)
        else:
            self.labels.append('')
        if profile is not None:
            self.profiles.append(profile.get_colors(self.pdb))
        else:
            self.profiles.append(self.pdb.get_colorby_sequence())

    def add_lines(self, view, i, j, color, viewer):
        pdb = self.pdb
        xi, yi, zi = pdb.get_xyz_coord(i)
        xj, yj, zj = pdb.get_xyz_coord(j)
        cylinder_specs = {"start": {"x": xi, "y": yi, "z": zi},
                          "end":  {"x": xj, "y": yj, "z": zj},
                          "radius": 0.5,
                          "fromCap": 2,
                          "toCap": 2,
                          "color": color}
        view.addCylinder(cylinder_specs, viewer=viewer)

    def plot_data(self, view, viewer, index):
        window = self.ij_windows[index]
        for i, j, color in zip(*self.ijs[index]):
            color = "0x"+mpc.rgb2hex(color)[1:]
            for w in range(window):
                io = i+w
                jo = j+window-1-w
                if io in self.pdb.validres and jo in self.pdb.validres:
                    self.add_lines(view, io, jo, color, viewer)

    def set_colors(self, view, viewer, index):
        colors = self.profiles[index]
        color_selector = {}
        for res in self.pdb.pdb.get_residues():
            res = res.get_id()
            color = colors[res[1]-1]
            if color in color_selector.keys():
                color_selector[color].append(res[1])
            else:
                color_selector[color] = [res[1]]
        for color in color_selector.keys():
            selector = {'resi': color_selector[color]}
            style = {"cartoon": {"color": color, "opacity": 0.8}}
            view.setStyle(selector, style, viewer=viewer)
        selector = {'resi': self.pdb.validres, 'invert': 'true'}
        style = {"cross": {"hidden": "true"}}
        view.setStyle(selector, style, viewer=viewer)

    def set_view(self, view):
        with open(self.pdb.path, 'r') as pdb_file:
            pdb_str = pdb_file.read()
        view.addModel(pdb_str, 'pdb')
        view.setStyle({'chain': 'A'}, {"cartoon": {'color': 'grey'}})
        view.zoomTo()
        return view

    def make_plot(self, view=None):
        rows, cols = get_rows_columns(len(self.ijs))
        if view is None:
            view = py3Dmol.view(viewergrid=(rows, cols),
                                width=400*rows, height=400*cols)
            view = self.set_view(view)
        if self.profiles is not None:
            colorby = "profile"
        else:
            colorby = "sequence"
        for i in range(len(self.ijs)):
            row = i // cols
            col = i % cols
            viewer = (row, col)
            self.set_colors(view, viewer, i)
            if self.ijs is not None:
                self.plot_data(view, viewer, i)
        return view
