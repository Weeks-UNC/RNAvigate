import py3Dmol
import matplotlib.colors as mpc
from .plots import Plot
import matplotlib.pyplot as plt


class Mol(Plot):
    def __init__(self, num_samples, pdb, width=400, height=400,
                 background_alpha=1):
        self.pdb = pdb
        self.length = num_samples
        self.rows, self.columns = self.get_rows_columns()
        view = py3Dmol.view(viewergrid=(self.rows, self.columns),
                            width=width*self.columns, height=height*self.rows)
        view.setBackgroundColor('white', background_alpha)
        with open(self.pdb.path, 'r') as pdb_file:
            pdb_str = pdb_file.read()
        view.addModel(pdb_str, 'pdb')
        view.setStyle({"cartoon": {'color': 'spectrum'}})
        view.zoomTo()
        self.view = view
        self.i = 0
        self.pass_through = ["nt_color", "atom", "title"]

    def get_figsize(self):
        pass

    def get_viewer(self, i=None):
        if i is None:
            i = self.i
        row = i // self.columns
        col = i % self.columns
        return (row, col)

    def plot_data(self, interactions, profile, label, nt_color="sequence",
                  atom="O2'", title=True):
        viewer = self.get_viewer()
        if interactions is not None:
            self.plot_interactions(viewer, interactions, atom)
            _, ax = plt.subplots(1, figsize=(6, 2))
            self.view_colormap(ax, interactions)
        self.set_colors(viewer, profile, nt_color)
        if title:
            self.view.addLabel(label,
                               {"position": {"x": 0, "y": 0, "z": 0},
                                "useScreen": True,
                                "fontColor": "black",
                                "backgroundOpacity": 0.0,
                                "fontSize": 28}, viewer=viewer)
        self.i += 1

    def add_lines(self, i, j, color, viewer, atom):
        pdb = self.pdb
        if i not in pdb.validres or j not in pdb.validres:
            return
        i += self.pdb.offset
        j += self.pdb.offset
        xi, yi, zi = pdb.get_xyz_coord(i, atom)
        xj, yj, zj = pdb.get_xyz_coord(j, atom)
        cylinder_specs = {"start": {"x": xi, "y": yi, "z": zi},
                          "end":  {"x": xj, "y": yj, "z": zj},
                          "radius": 0.5,
                          "fromCap": 2,
                          "toCap": 2,
                          "color": color}
        self.view.addCylinder(cylinder_specs, viewer=viewer)

    def plot_interactions(self, viewer, interactions, atom):
        window = interactions.window
        for i, j, color in zip(*interactions.get_ij_colors()):
            color = "0x"+mpc.rgb2hex(color)[1:]
            for w in range(window):
                io = i+w
                jo = j+window-1-w
                self.add_lines(io, jo, color, viewer, atom)

    def set_colors(self, viewer, profile, nt_color):
        colors = self.pdb.get_colors(nt_color, profile=profile)
        color_selector = {}
        valid_pdbres = []
        for i, res in enumerate(self.pdb.validres+self.pdb.offset):
            res = int(res)
            valid_pdbres.append(res)
            color = colors[i]
            if color in color_selector.keys():
                color_selector[color].append(res)
            else:
                color_selector[color] = [res]
        for color in color_selector.keys():
            selector = {'chain': self.pdb.chain, 'resi': color_selector[color]}
            style = {"cartoon": {"color": color}}
            self.view.setStyle(selector, style, viewer=viewer)
        selector = {'chain': self.pdb.chain,
                    'resi': valid_pdbres, 'invert': 'true'}
        style = {"cross": {"hidden": "true"}}
        self.view.setStyle(selector, style, viewer=viewer)

    def hide_cylinders(self):
        resns = ['A', 'G', 'C', 'U']
        atoms = ['N1', 'N1', 'N3', 'N3']
        for resn, atom in zip(resns, atoms):
            self.view.setStyle(
                {'resn': resn, 'atom': atom},
                {'cross': {'hidden': 'true'}})

    def png(self):
        '''output png image of viewer, which must already be instantiated'''
        script = '''<script>
            var pngdata = viewer_{0}.pngURI()
            </script>'''.format(self.view.uniqueid)
        print("To view png in notebook, type:\n"
              "IPython.display.publish_display_data({'text/html': plot.png()})\n"
              "Then, to save: right click image and click 'save as'\n"
              "Currently not working correctly in VSCode, only Jupyter")
        return script
