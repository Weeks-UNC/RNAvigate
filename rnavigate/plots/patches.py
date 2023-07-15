# matplotlib patches defined by data units

import matplotlib as mpl

class LineDataUnits(mpl.lines.Line2D):
    """Line2D subclass. The only difference is that linewidth is computed
    in data coordinates and updated any time the figure size changes.
    """
    def __init__(self, *args, **kwargs):
        _lw_data = kwargs.pop("linewidth", 1) 
        super().__init__(*args, **kwargs)
        self._lw_data = _lw_data

    def _get_lw(self):
        if self.axes is not None:
            ppd = 72./self.axes.figure.dpi
            trans = self.axes.transData.transform
            return ((trans((1, self._lw_data))-trans((0, 0)))*ppd)[1]
        else:
            return 1

    def _set_lw(self, lw):
        self._lw_data = lw

    _linewidth = property(_get_lw, _set_lw)

def plot_circle_markers(ax, centers, radii, colors):
    circles = []
    for center, radius, color in zip(centers, radii, colors):
        circles.append(mpl.patches.Circle(center, radius=radius, linewidth=0))
    c = mpl.collections.PatchCollection(circles, match_original=True)
    ax.add_collection(c)
