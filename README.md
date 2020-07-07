Jupyter Notebook Plotting Tools
==============================================================================
Python plotting functions for making pretty figures for common lab experiments.
If you'd like me to add any functions to this, let me know by openning up an
issue in the issues tab.

In this README, I will simply list out the functions associated with each of
the modules contained in this package. For more in-depth exploration of the
functions and example usage, see the
[wiki](https://github.com/Weeks-UNC/JNBTools/wiki).

plottingTools.py
------------------------------------------------------------------------------
```python
readHistograms(logfile):
plotSkyline(axis, profile, label=None, column='Reactivity_profile'):
addSeqBar(axis, profile, yvalue=-0.017):
getWidth(sample):
```


shapemapperPlots.py
------------------------------------------------------------------------------
```python
plotProfile(axis, sample, name)
plotMutationRates(axis, sample)
plotDepth(axis, sample)
```


JNBarcPlot.py
------------------------------------------------------------------------------
```python
plotIgnoredCorrs(ax, title, logfile)
arcPlot(ct=False, fasta=False, refct=False, probability=False, ringz=False,
        ringsig=False, pairmap=False, compare_pairmap=False, ntshape=False,
        dmsprofile=False, bottom=False, title=False, showGrid=False,
        bound=False, filternc=False, profile=False):
```


BMTools
------------------------------------------------------------------------------
```python
class BM():
    def __init__(self, reactivities):
    def plotBMmodel(self, axis):
```
