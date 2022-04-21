import py3Dmol
import matplotlib.colors as mpc
from .plots import Plot
import matplotlib.pyplot as plt


class Mol(Plot):
    def __init__(self, num_samples, pdb):
        self.pdb = pdb
        self.length = num_samples
        self.rows, self.columns = self.get_rows_columns()
        view = py3Dmol.view(viewergrid=(self.rows, self.columns),
                            width=800*self.columns, height=800*self.rows)
        with open(self.pdb.path, 'r') as pdb_file:
            pdb_str = pdb_file.read()
        view.addModel(pdb_str, 'pdb')
        view.setStyle({"cartoon": {'color': 'spectrum', 'opacity': 0.8}})
        view.zoomTo()
        self.view = view
        self.i = 0
        self.pass_through = ["cmap", "nt_color"]

    def get_figsize(self):
        pass

    def get_viewer(self, i=None):
        if i is None:
            i = self.i
        row = i // self.columns
        col = i % self.columns
        return (row, col)

    def plot_data(self, ij, profile, label, cmap=None, nt_color="sequence"):
        viewer = self.get_viewer()
        if ij is not None:
            self.plot_ij(viewer, ij, cmap=cmap)
            _, ax = plt.subplots(1, figsize=(6, 2))
            self.view_colormap(ax, ij, cmap=cmap)
        self.set_colors(viewer, profile, nt_color)
        self.view.addLabel(label,
                           {"position": {"x": 400, "y": 50, "z": 0},
                            "useScreen": True, "alignment": "center",
                            "fontColor": "black", "backgroundColor": "white",
                            "fontSize": 28}, viewer=viewer)
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

    def plot_ij(self, viewer, ij, cmap=None):
        window = ij.window
        for i, j, color in zip(*ij.get_ij_colors(cmap=cmap)):
            color = "0x"+mpc.rgb2hex(color)[1:]
            for w in range(window):
                io = i+w
                jo = j+window-1-w
                self.add_lines(io, jo, color, viewer)

    def set_colors(self, viewer, profile, nt_color):
        if nt_color in ["sequence", "profile"]:
            colors = self.pdb.get_colors(nt_color, profile=profile)
        elif nt_color == "position":
            return
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
