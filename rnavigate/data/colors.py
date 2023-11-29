import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mpc


class ScalarMappable(cm.ScalarMappable):
    def __init__(self, cmap, normalization, values, title="", tick_labels=None,
                 **cbar_args):
        cmap = self.get_cmap(cmap)
        norm = self.get_norm(normalization, values, cmap)
        super().__init__(norm, cmap)
        self.rnav_norm = normalization
        self.rnav_vals = values
        self.rnav_cmap = cmap(np.arange(cmap.N))
        self.cbar_args = cbar_args
        self.tick_labels = tick_labels
        self.title = title

    def is_equivalent_to(self, cmap2):
        return ((self.rnav_norm == cmap2.rnav_norm)
                and (self.rnav_vals == cmap2.rnav_vals)
                and (len(self.rnav_cmap) == len(cmap2.rnav_cmap))
                and mpc.same_color(self.rnav_cmap, cmap2.rnav_cmap)
                and (self.tick_labels == cmap2.tick_labels)
                and (self.cbar_args == cmap2.cbar_args)
                and (self.title == cmap2.title))

    def values_to_hexcolors(self, values, alpha=1.0):
        colors = super().to_rgba(x=values, alpha=alpha)
        return np.array([mpc.to_hex(c, keep_alpha=True) for c in colors])

    def get_norm(self, normalization, values, cmap):
        if normalization == 'min_max':
            return mpc.Normalize(values[0], values[1])
        elif normalization == '0_1':
            return mpc.Normalize()
        elif normalization == 'none':
            return mpc.NoNorm()
        elif normalization == 'bins':
            if len(values) == (cmap.N - 1):
                return mpc.BoundaryNorm(values, cmap.N, extend='both')
            if len(values) == (cmap.N + 1):
                return mpc.BoundaryNorm(values, cmap.N)

    def get_cmap(self, cmap):
        """Given a matplotlib color, list of colors, or colormap name, return
        a colormap object

        Args:
            cmap (str | tuple | float | list): A valid mpl color, list of valid
                colors or a valid colormap name

        Returns:
            matplotlib colormap: listed colormap matching the input
        """
        if mpc.is_color_like(cmap):
            cmap = mpc.ListedColormap([cmap])
            cmap.set_bad('grey')
            return cmap
        elif (isinstance(cmap, list)
                and all(mpc.is_color_like(c) for c in cmap)):
            cmap = mpc.ListedColormap(cmap)
            cmap.set_bad('grey')
            return cmap
        try:
            return cm.get_cmap(cmap)
        except ValueError as exception:
            print("cmap must be one of: valid mpl color, list of mpl colors, or "
                f"mpl colormap:\n{str(cmap)}")
            raise exception
