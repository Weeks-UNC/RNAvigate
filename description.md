# Class and methods descriptions
## ./starmapper.py - import starmapper as MaP

<details><summary>Class: Sample()</summary>
<details><summary>Initialization</summary>

  - __init__()
    - For each parameter that gets passed during Sample initialization, a data
      object is created and added to the "data" dictionary property.
      - e.g. Sample(ct="ctfile.ct") -> Sample.data["ct"] = CT(file="ctfile.ct")
    - If dance_prefix is provided, dance components are stored as a list of
        Sample objects in the property "Sample.dance". This will include PAIR,
        RInG, CT, and bp probabilities if files are found matching the prefix.
    - Parameters:
      - sample = a name or label given to this sample for figure legends
      - fasta = a reference sequence, required for JuMP data
      - profile = a shapemapper profile.txt output file
      - ct = a ct structure to be associated with this data
      - compct = a second ct structure
      - ss = a secondary structure file (.xrna, .varna, .cte, or .nsd)
      - log = a shapemapper log output file (_log.txt)
      - rings = a correlation output file from ringmapper
      - deletions = a deletions count output file from SHAPE-JuMP
      - pairs = a PAIRs output file from pairmapper
      - pdb = a pdb structure, must be standard pdb format
      - pdb_kwargs = dictionary containing arguments for parsing pdb files
        - chain = name of the chain of interest, should be a single character
        - fasta = a reference sequence, required only if header is missing
        - offset = numeric position for 1st nt, required only if no header
      - probs = a base-pairing probability dotplot file from ProbabilityPlot
      - dance_prefix = the file name prefix for all dance-mapper output files
</details>
<details><summary>Internal methods</summary>

  - init_dance(self, prefix)
    - Called during Sample initialization. This function looks for files
      matching the given dance file prefix and creates a list of sample objects
      containing data about each component of the dance model.
  - get_data(self, key)
    - looks for key (string) to return the corresponding property.
  - get_data_list(self, *keys)
    - calls get_data on each key and packs result into a list.
  - filter_ij(self, ij, fit_to, **kwargs)
    - convenience function for ij.filter(), Sample.filter_ij('rings', 'ct') is
      equivalent to Sample.data["rings"].filter(fit_to=Sample.data["ct"])
  - dance_filter(self, filterneg=True, cdfilter=15, sigfilter=23, ssfilter=True)
    - applies filters appropriate for plotting dance rings and pairs together.
</details>
<details><summary>Plotting methods</summary>

  - make_*plot() - makes a single plot from Sample data
  - make_*plot_multifilter() - makes multiple plots from sets of filters
    - Replace *plot with skyline, shapemapper, ap, ss, mol, heatmap, circle or qc
</details>
</details>

<details><summary>Functions</summary>
<details><summary>Function: create_code_button()</summary>

  When used within a Jupyter notebook, this will embed an HTML button which 
  toggles hiding/showing code blocks. This is especially useful for creating a
  report of your analysis exported to HTML.
</details>
<details><summary>Function: array_*plot()</summary>

  Replace *plot with qc, skyline, ap, ss, mol, heatmap, circle, or linreg. These
  take a list of samples and makes one plot per sample. Makes it easy to compare
  similar data between samples.
  
  Sample.make_plot() == MaP.array_plot([Sample])
</details>
</details>

## ./plots/
<details><summary>plots.py - Abstract class: Plot()</summary>

This is an abstract class that defines some generally useful methods and what
properties a plot object should have.
### Methods:
- __init__() from the number of samples, makes a figure and axis grid.
- get_ax() return the axis corresponding to the *n*th sample
- add_sample() retrieves data and passes it on to plot_data()
- view_colormap() creates an appropiate colormap for given ij data.
- get_rows_columns() returns # of rows an columns for a given # of samples
- add_sequence() adds a sequence bar along the x-axis of the given plot
- abstract get_figsize(): each figure must know how big to be.
- abstract plot_data(): each figure must know how to plot the given data
</details>
<details><summary>arc.py - Class: AP()</summary>

This creates a grid of arcPlots.
### Methods
- plot_data() plots the given data on the current axis, then moves to next axis.
- add_patches() draws arcs as a patch_collection.
- add_title()
- get_figsize()
- plot_profile() draws the mid-plot reactivity bar chart
</details>
<details><summary>circle.py - Class: Circle()</summary>

Creates a grid of circle plots, similar to an arcPlot, but arranged in a circle.
### Methods
- plot_data() plots the data on the current axis, then moves to the next.
- get_figsize()
- add_patches() draws arcs as a patch collection
</details>
<details><summary>heatmap.py - Class: Heatmap()</summary>

--
A nt x nt grid that plots each ij data point as a pixel. Useful for viewing very
dense data.
### Methods
- plot_data() adds a structure as contour map and ij data as heatmap.
- get_figsize()
- plot_contour_distances() draws the contour plot
- plot_heatmap_data() draws the heatmap
</details>
<details><summary>Linear Regression</summary>

---

### linreg.py - Class LinReg()
A classic grid plot. Showing separation of paired and unpaired reactivity along
diagonal and linear regression between each passed sample in the other positions.
### Methods
- get_rows_columns()
- get_figsize()
- plot_data()
- plot_regression()
- plot_kde()

---

</details>
<details><summary>3-D Molecules</summary>

---

### mol.py - Class Mol()
An interactive 3D RNA structure colored by reactivity, with ij data plotted as
cylinders. Built with py3DMol (3DMol.js).
### Methods
- get_figsize()
- get_viewer() similar to get_ax() but for a 3DMol view object
- plot_data()
- add_lines()
- plot_ij()
- set_colors()

---

</details>
<details><summary>Quality Control</summary>

---

### qc.py - Class: QC()
Quality control metrics including distributions of mutations per molecule, read
lengths, and reactivities.
### Methods
- 
</details>
<details><summary>skyline.py - Class: Skyline()</summary></details>
<details><summary>sm.py - Class: SM()</summary></details>
<details><summary>ss.py - Class: SS()</summary></details>

## ./data/
<details><summary>data.py - Abstract class: Data()</summary></details>
<details><summary>ct.py - Class: CT()</summary></details>
<details><summary>dp.py - Class: DP()</summary></details>
<details><summary>ij.py - Class: IJ()</summary></details>
<details><summary>log.py - Class: Log()</summary></details>
<details><summary>pdb.py - Class: PDB()</summary></details>
<details><summary>profile.py - Class: Profile()</summary></details>
