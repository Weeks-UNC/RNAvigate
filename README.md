Jupyter Notebook Plotting Tools
==============================================================================
Python plotting functions for making pretty figures for common lab experiments.
If you'd like me to add any functions to this, let me know by openning up an
issue in the issues tab.

Templates and Code Examples
------------------------------------------------------------------------------
Take a look at the ipynb-template for example plots and the code to make them.

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
