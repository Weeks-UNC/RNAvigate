# Class and methods descriptions
## ./plotmapper.py - import plotmapper as MaP

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
<details><summary>Function: create_code_button()</summary>

</details>
</details>

## ./plots
<details><summary>plots.py - Abstract class: Plot()</summary></details>
<details><summary>arc.py - Class: AP()</summary></details>
<details><summary>circle.py - Class: Circle()</summary></details>
<details><summary>heatmap.py - Class: Heatmap()</summary></details>
<details><summary>linreg.py - Class: LinReg()</summary></details>
<details><summary>mol.py - Class: Mol()</summary></details>
<details><summary>qc.py - Class: QC()</summary></details>
<details><summary>skyline.py - Class: Skyline()</summary></details>
<details><summary>sm.py - Class: SM()</summary></details>
<details><summary>ss.py - Class: SS()</summary></details>

## ./data
<details><summary>data.py - Abstract class: Data()</summary></details>
<details><summary>ct.py - Class: CT()</summary></details>
<details><summary>dp.py - Class: DP()</summary></details>
<details><summary>ij.py - Class: IJ()</summary></details>
<details><summary>log.py - Class: Log()</summary></details>
<details><summary>pdb.py - Class: PDB()</summary></details>
<details><summary>profile.py - Class: Profile()</summary></details>
