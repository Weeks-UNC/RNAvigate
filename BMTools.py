import plottingTools as pt
import pandas as pd


class BM():
    """A class for handling reactivities.txt files"""
    def __init__(self, reactivities):
        with open(reactivities) as inf:
            header1 = inf.readline().strip().split()
            header2 = inf.readline().strip().split()
        self.components = int(header1[0])
        self.p = header2[1:]
        colnames = ["Nucleotide", "Sequence"]
        for i in range(self.components):
            colnames.append("nReact"+str(i))
            colnames.append("Raw"+str(i))
            colnames.append("blank"+str(i))
        colnames.append("Background")
        self.reactivities = pd.read_csv(reactivities, sep='\t', header=2,
                                        names=colnames)

    def plotBMmodel(self, axis):
        """Plot all profiles as skyline. Legend includes population values."""
        for i in range(self.components):
            pt.plotSkyline(axis, self.reactivities,
                           label="{}: {}".format(i, self.p[i]),
                           column="nReact{}".format(i))
        pt.addSeqBar(axis, self.reactivities, yvalue=-0.1)
        axis.legend(title="Component: Population")
