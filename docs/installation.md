Installation
============

- [Anaconda installation](#anaconda-installation)
- [Manual installation](#manual-installation)
- [UNC Longleaf installation](#unc-longleaf-installation)

Anaconda installation
---------------------

RNAvigate is not available via pip or conda. This is coming soon. However, I
have included a conda environment file. If you are using conda, follow this
guide for the easiest installation.

1. Open a terminal and navigate to the directory you'll be saving RNAvigate.
2. Clone the github repository.
3. Navigate into the RNAvigate directory.
4. Create and activate a conda enviroment using the included env.yaml file
5. Install the RNAvigate package in the new conda environment.

```bash
git clone https://github.com/Weeks-UNC/RNAvigate.git
cd RNAvigate
conda env create -f env.yaml
conda activate rnavigate
conda develop .
```

If you are planning to use Jupyter Notebooks to do your analyses (recommended),
you will need to make the new conda environment available to Jupyter.

```bash
python -m ipykernel install --user --name=rnavigate
```

To test installation, open and run one of the example notebooks.

Manual installation
-------------------

This section is bit vague for now. I will make RNAvigate available as a pip
package soon, and this section will hopefully be obsolete.

To install manually, you will need the following dependencies in your python
environment:
- scipy >= 1.6.2
- py3dmol >= 0.8.0
- matplotlib >= 3.3.4
- python >= 3.9.2
- seaborn >= 0.11.1
- biopython >= 1.78
- numpy >= 1.20.1
- pandas >= 1.2.4
- scikit-learn >= 1.0.2

Download the package from [Github](https://github.com/Weeks-UNC/RNAvigate). Add
this directory your $PYTHONPATH, and you are good to go.

UNC Longleaf installation
-------------------------

Log into Longleaf using ssh and run the following code from the command line.

```bash
cd $HOME
git clone https://github.com/Weeks-UNC/RNAvigate.git
```

Next we need to make sure that you have an Anaconda environment that includes
all dependencies and that that environment is available to Jupyter. This does
not change your default modules. They will be restored next time you log in.

```bash
module rm python pymol pyrosetta
module load anaconda/2019.10
conda env create -f /proj/kweeks/bin/RNAvigate/env.yaml
conda activate rnavigate
conda develop /proj/kweeks/bin/RNAvigate/
python -m ipykernel install --user --name=rnavigate
```

If this occured without errors, exit longleaf and go to
[UNC's OpenOnDemand Service](https://ondemand.rc.unc.edu/). OpenOnDemand is
UNC's platform for using interactive programs within Longleaf.

1. Log in using your ONYEN, and start a jupyter notebook.

   1. Click on interactive apps, under servers, click on Jupyter Notebook.
   2. Enter the number of hours you will be working, 1 CPU, other fields blank.
   3. DONT FORGET TO SAVE BEFORE YOUR NOTEBOOK SHUTS DOWN.

2. Navigate to your data directory, and open a new notebook using the
  "rnavigate" option.
4. Type `import rnavigate as MaP` into the first code cell and run it.
5. Your ready to starting exploring your data!
