import py3Dmol
import matplotlib.colors as mpc
from .plots import Plot
import matplotlib.pyplot as plt


class Mol(Plot):
    def __init__(self, num_samples, pdb):
        self.pdb = pdb
        self.rows, self.columns = self.get_rows_columns(num_samples)
        view = py3Dmol.view(viewergrid=(self.rows, self.columns),
                            width=800*self.columns, height=800*self.rows)
        with open(self.pdb.path, 'r') as pdb_file:
            pdb_str = pdb_file.read()
        view.addModel(pdb_str, 'pdb')
        view.setStyle({"cartoon": {'color': 'grey'}})
        view.zoomTo()
        self.view = view
        self.i = 0
        self.pass_through = []

    def get_figsize(self):
        pass

    def get_viewer(self, i=None):
        if i is None:
            i = self.i
        row = i // self.columns
        col = i % self.columns
        return (row, col)

    def plot_data(self, ij, profile, label):
        viewer = self.get_viewer()
        if ij is not None:
            self.plot_ij(viewer, ij)
            _, ax = plt.subplots(1, figsize=(6, 2))
            self.view_colormap(ax, ij)
        self.set_colors(viewer, profile)
        print(f"viewer: {viewer}, {label}")
        self.i += 1

    def add_lines(self, i, j, color, viewer):
        pdb = self.pdb
        if i not in pdb.validres or j not in pdb.validres:
            return
        i += self.pdb.offset
        j += self.pdb.offset
        xi, yi, zi = pdb.get_xyz_coord(i)
        xj, yj, zj = pdb.get_xyz_coord(j)
        cylinder_specs = {"start": {"x": xi, "y": yi, "z": zi},
                          "end":  {"x": xj, "y": yj, "z": zj},
                          "radius": 0.5,
                          "fromCap": 2,
                          "toCap": 2,
                          "color": color}
        self.view.addCylinder(cylinder_specs, viewer=viewer)

    def plot_ij(self, viewer, ij):
        window = ij.window
        for i, j, color in zip(*ij.get_ij_colors()):
            color = "0x"+mpc.rgb2hex(color)[1:]
            for w in range(window):
                io = i+w
                jo = j+window-1-w
                self.add_lines(io, jo, color, viewer)

    def set_colors(self, viewer, profile):
        if profile is None:
            colors = self.pdb.get_colorby_sequence()
        else:
            colors = profile.get_colors(self.pdb)
        color_selector = {}
        valid_pdbres = []
        for res in self.pdb.validres:
            res_off = res + self.pdb.offset
            valid_pdbres.append(res_off)
            color = colors[res-1]
            if color in color_selector.keys():
                color_selector[color].append(res_off)
            else:
                color_selector[color] = [res_off]
        for color in color_selector.keys():
            selector = {'chain': self.pdb.chain, 'resi': color_selector[color]}
            style = {"cartoon": {"color": color, "opacity": 0.8}}
            self.view.setStyle(selector, style, viewer=viewer)
        selector = {'chain': self.pdb.chain,
                    'resi': valid_pdbres, 'invert': 'true'}
        style = {"cross": {"hidden": "true"}}
        self.view.setStyle(selector, style, viewer=viewer)
