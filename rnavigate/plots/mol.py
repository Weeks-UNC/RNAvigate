import py3Dmol
import matplotlib.colors as mpc
from .plots import Plot
import matplotlib.pyplot as plt


class Mol(Plot):
    def __init__(self, num_samples, pdb, width=400, height=400,
                 background_alpha=1, rotation=None, orientation=None,
                 rows=None, cols=None):
        self.pdb = pdb
        self.length = num_samples
        self.rows, self.columns = self.get_rows_columns(rows=rows, cols=cols)
        view = py3Dmol.view(viewergrid=(self.rows, self.columns),
                            width=width*self.columns, height=height*self.rows)
        view.setBackgroundColor('white', background_alpha)
        with open(self.pdb.path, 'r') as pdb_file:
            pdb_str = pdb_file.read()
        if self.pdb.path.split('.')[-1] == "pdb":
            view.addModel(pdb_str, 'pdb')
        elif self.pdb.path.split('.')[-1] == "cif":
            view.addModel(pdb_str, 'mmcif')
        view.setStyle({"cartoon": {'color': 'spectrum'}})
        view.zoomTo({"chain": self.pdb.chain})
        if orientation is not None:
            view.setView(orientation)
        elif rotation is not None:
            for key in rotation:
                view.rotate(rotation[key], key)
        self.view = view
        self.i = 0
        self.pass_through = [
            "nt_color",
            "atom",
            "title",
            "get_orientation",
            "colorbar"]

    def get_figsize(self):
        pass

    def get_viewer(self, i=None):
        if i is None:
            i = self.i
        row = i // self.columns
        col = i % self.columns
        return (row, col)

    def plot_data(self, interactions, profile, label, nt_color="grey",
                  atom="O2'", title=True, get_orientation=False,
                  colorbar=True, viewer=None):
        if viewer is None:
            viewer = self.get_viewer()
        if get_orientation:
            self.get_orientation()
        if interactions is not None:
            self.plot_interactions(viewer, interactions, atom)
            if colorbar:
                self.view_colormap(interactions=interactions)
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
        if not pdb.is_valid_idx(seq_idx=i) or not pdb.is_valid_idx(seq_idx=j):
            return
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
        for res in self.pdb.pdb_idx:
            res = int(res)
            valid_pdbres.append(res)
            color = colors[self.pdb.get_seq_idx(res)-1]
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

    def save(self):
        '''output png image of viewer, which must already be instantiated'''
        print("To save, orient the interactive plot to the view you'd like\n"
              "save, then run plot.save() in a new cell. The resulting image\n"
              "will be saveable as a png file")
        return self.view.png()

    def get_orientation(self):
        self.view.setClickable(
            {},
            'true',
            '''function(atom,viewer,event,container){
                viewer.addLabel(viewer.getView().map(function(each_element){return each_element.toFixed(2)}),
                {
                    position: {x:0,y:0,z:0},
                    useScreen: true,
                    backgroundColor: 'black'
                });
            }'''
        )
