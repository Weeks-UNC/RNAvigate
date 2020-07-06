Jupyter Notebook Plotting Tools
==============================================================================
Python plotting functions for making pretty figures for common lab experiments.
If you'd like me to add any functions to this, let me know by openning up an
issue in the issues tab.

Table of Contents
------------------------------------------------------------------------------

Hiding code blocks
------------------------------------------------------------------------------
You can put this code in the first code block of your notebook. This will
hide your code by default and produce a button which will toggle the code
blocks off and on. This is useful for making a report more readable, and having
the option to show the code if it is needed.

```
from IPython.display import HTML

HTML('''<script>
code_show=true;
function code_toggle() {
 if (code_show){
 $('div.input').hide();
 } else {
 $('div.input').show();
 }
 code_show = !code_show
}
$( document ).ready(code_toggle);
</script>
<form action="javascript:code_toggle()"><input type="submit" value="Click here to toggle on/off the raw code."></form>''')
```


Loading in modules
------------------------------------------------------------------------------
Modules in your path can be imported easily. You will also need matplotlib,
seaborn, and pandas. Import only the modules that you'll be using.

```
import plottingTools as pt
import JNBarcPlot as ap
import BMTools as bm
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
```

You will also want to set up your plotting style defaults at the top of your
notebook. Here are the ones that I use:

```
sns.set_style("ticks")
sns.set_context("talk")
colors = ['#a100ffff', '#edc600ff', '#0092edff', '#ff8300ff', '#ff48e9ff', '#3fd125ff']
sns.set_palette(colors)
```

If you are working from the longleaf jupyter hub, you will first need to add
JNBTools, arcPlot, and RNATools to your path. Edit this code to point to the
right directory. These should go before the rest of your import statements.

```
import sys
sys.path.append("/nas/longleaf/home/<ONYEN>/JNBTools")
sys.path.append("/nas/longleaf/home/<ONYEN>/arcPlot")
sys.path.append("/nas/longleaf/home/<ONYEN>/RNATools")
```


Loading in data
------------------------------------------------------------------------------
These functions all work using pandas DataFrames. Therefor, loading a file is an
easy one-liner.
```
# for shapemapper profiles
shapedata = pd.read_csv('example_RNA_profile.txt', sep='\t')

# for pairmapper/ringmapper correlation files
ringdata = pd.read_csv('example_RNA-pairmap.txt', sep='\t', header=1)

# for BM profiles, these have their own class and methods.
bmdata = bm.BM('example_RNA_profiles.txt')
```


Example plot: Tried-and-true shapemapper plots
------------------------------------------------------------------------------
Here is an example of how to make the standard shapemapper2 plots that we are
all used to seeing.
```
fig, ax = plt.subplots(3, 1, figsize=(pt.getWidth(sample), 15))
pt.plotProfile(ax[0], shapedata, 'Sample Name')
pt.plotMutationRates(ax[1], shapedata)
pt.plotDepth(ax[2], shapedata)
```

![Standard ShapeMapper Plots](./example-plots/standard-shapemapper-plots.png)


Notes About Setup
------------------------------------------------------------------------------
Running this code locally requires that you have ArcPlot and RNATools2 in the
parent directory. I'm writing this code with Python 3 notebooks in mind.
However, you may find that it is backwords compatible with Python 2. If that
isn't the case, then you will need to convert ArcPlot and RNATools2 to be
compatible with Python 3. To do this, navigate to the directories that hold
these scripts and run the following code.

```
cd arcPlot
2to3 -w arcPlot.py
2to3 -w pmanalysis.py
cd ../RNATools
2to3 -w RNAtools2.py
```

Use in Longleaf
------------------------------------------------------------------------------
You can access a jupyter hub to create notebooks within longleaf by going
[here](https://longleaf-jupyter.its.unc.edu/hub/login).
