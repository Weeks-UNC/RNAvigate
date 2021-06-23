from plots import *
import py3Dmol
import matplotlib.colors as mpc


class Mol():
    def __init__(self, pdbs, profiles=None, ijs=None, labels=None):
        self.pdbs = pdbs
        self.labels = labels
        self.profiles = profiles
        self.ijs = ijs

    def add_lines(self, view, i, j, color, viewer, index):
        pdb = self.pdbs[index]
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
        pdb = self.pdbs[index]
        ij = self.ijs[index]
        ij_colors = ij.get_ij_colors()
        for i, j, color in zip(*ij_colors):
            color = mpc.rgb2hex(color)
            color = "0x"+color[1:]
            window = ij.window
            for w in range(window):
                io = i+w
                jo = j+window-1-w
                if io in pdb.validres and jo in pdb.validres:
                    self.add_lines(view, io, jo, color, viewer, index)

    def set_colors(self, view, colorby, viewer, index):
        pdb = self.pdbs[index]
        if colorby == "profile":
            colors = self.profiles[index].get_colors(pdb)
        if colorby == "sequence":
            colors = pdb.get_colorby_sequence()
        color_selector = {}
        for res in pdb.pdb.get_residues():
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
        selector = {'resi': pdb.validres, 'invert': 'true'}
        style = {"cross": {"hidden": "true"}}
        view.setStyle(selector, style, viewer=viewer)

    def set_view(self, view):
        with open(self.pdbs[0].path, 'r') as pdb_file:
            pdb_str = pdb_file.read()
        view.addModel(pdb_str, 'pdb')
        view.setStyle({'chain': 'A'}, {"cartoon": {'color': 'grey'}})
        view.zoomTo()
        return view

    def make_plot(self, view=None):
        rows, cols = get_rows_columns(len(self.pdbs))
        if view is None:
            view = py3Dmol.view(viewergrid=(rows, cols),
                                width=400*rows, height=400*cols)
            view = self.set_view(view)
        if self.profiles is not None:
            colorby = "profile"
        else:
            colorby = "sequence"
        for i in range(len(self.pdbs)):
            row = i // cols
            col = i % cols
            viewer = (row, col)
            self.set_colors(view, colorby, viewer, i)
        if self.ijs is not None:
            self.plot_data(view, viewer, i)
        return view
