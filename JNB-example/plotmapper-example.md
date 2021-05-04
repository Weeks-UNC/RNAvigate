
Example Notebook for using plotmapper
=====================================
plotmapper.py will replace and combine the functionality of JNBarcPlot.py, plottingTools.py, secondaryStructure.py, and shapemapperPlots.py using a unified class structure (Sample for a single MaP sample, Experiment for multiple samples). This will allow much higher level figure creation, but it is not quite yet ready.

Below are examples of the high level functions available for plotting a variety of MaP data.

Contents:
- [Notebook set-up](#notebook-set-up)
- [Initializing MaP sample](#initializing-map-sample)
- [High-level plotting functions](#high-level-plotting-functions)
- [ShapeMapper QC](#shapemapper-qc)
- [Classic ShapeMapper Plots](#classic-shapemapper-plots)
- [Skyline plots](#skyline-plots)
- [DANCE-MaP reactivity skyline](#dance-map-reactivity-skyline)
- [Arc Plots](#arc-plots)
- [Secondary Structure](#secondary-structure)

Notebook set-up
---------------


```python
# This sets plots to display in-line by default
%matplotlib inline

# Import module, for high-level functions, no additional modules are needed
import plotmapper as MaP

# Creates an HTML button that hides/shows code cells
# Useful for lab notebook reports and research updates
MaP.create_code_button()
```


<script>
                 code_show=true;
                 function code_toggle() {
                 if (code_show) {$('div.input').hide();}
                 else {$('div.input').show();}
                 code_show = !code_show
                 }
                 $( document ).ready(code_toggle);
                 </script>
                 <form action="javascript:code_toggle()">
                 <input type="submit" value="Hide/show raw code.">
                 </form>


Initializing MaP sample
-----------------------
If you have consistently named files, (which you should), you can use a function to create a dictionary of keyword arguments (kwargs). Then, "unpack" the dictionary using the double asterisk.


```python
path = 'data/'
def kwargs(sample):
    kwargs = {}
    kwargs["profile"] = path+sample+"_rnasep_profile.txt"
    kwargs["ctfile"] = path+"RNaseP.ct"
    kwargs["ss"] = path+"RC_CRYSTAL_STRUCTURE.xrna"
    kwargs["ringfile"] = path+sample+"-rnasep.corrs"
    kwargs["pairfile"] = path+sample+"-rnasep-pairmap.txt"
    kwargs["logfile"] = path+sample+"_shapemapper_log.txt"
    kwargs["dance_reactivities"] = path+sample+"_rnasep-reactivities.txt"
    return kwargs

example = MaP.Sample(sample="example", **kwargs("example"))
# equivalent to:
# example = MaP.Sample(sample="example",
#                      profile='data/example_rnasep_profile.txt',
#                      ctfile='data/RNaseP.ct',
#                      ss='data/RC_CRYSTAL_STRUCTURE.xrna',
#                      ringfile='data/example-rnasep.corrs',
#                      pairfile='data/exampleexample-rnasep-pairmap.txt',
#                      logfile='data/example_shapemapper_log.txt',
#                      dance_reactivities='data/example_rnasep-reactivities.txt')
```

    Note: T nucleotides have been recoded as U
    

High-level plottting functions
------------------------------

ShapeMapper QC
--------------


```python
example.make_log_qc()

```


![svg](plotmapper-example_files/plotmapper-example_7_0.svg)


Classic ShapeMapper Plots
-------------------------


```python
example.make_shapemapper()

```


![svg](plotmapper-example_files/plotmapper-example_9_0.svg)


Skyline Plots
-------------


```python
example.make_skyline()
# def make_skyline(self, column="Reactivity_profile"):
#     fig, ax = plt.subplots(1, figsize=self.get_skyline_figsize(1, 1))
#     self.plot_skyline(ax)
#     self.plot_sequence(ax)
#     ax.set(title="Raw Reactivity Profile",
#            xlim=[0, self.length],
#            xticks=range(0, self.length, 20))
#     ax.legend(title="Samples")
```


![svg](plotmapper-example_files/plotmapper-example_11_0.svg)


DANCE-MaP reactivities skyline
------------------------------


```python
example.make_dance_skyline()
```


![svg](plotmapper-example_files/plotmapper-example_13_0.svg)


Arc Plots
---------


```python
example.make_ap(type="pairs")
```


![svg](plotmapper-example_files/plotmapper-example_15_0.svg)


Secondary Structure
-------------------


```python
example.make_ss()
```


![svg](plotmapper-example_files/plotmapper-example_17_0.svg)

