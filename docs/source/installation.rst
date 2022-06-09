============
Installation
============

.. note::

   If you are using UNC's Longleaf computing cluster, skip ahead to "Longleaf
   installation".

1. Open a terminal and navigate to the directory you'll be saving RNAvigate.
2. Clone the github repository.
3. Navigate into the RNAvigate directory.
4. Create and activate a conda enviroment using the included env.yaml file
5. Install the RNAvigate package in the new conda environment.

.. code-block::
   :caption: Downloading the github repository

   git clone https://github.com/Weeks-UNC/RNAvigate.git
   cd RNAvigate
   conda env create -f env.yaml -n rnavigate
   source activate rnavigate
   pip install .

If you are planning to use Jupyter Notebooks to do your analyses (highly
recommended), you will need to make the new conda environment available to Jupyter.
.. code-block::
   :caption: Creating a Jupyter Notebook conda environment

   python -m ipykernel install --user --name=rnavigate

If this occured without errors, your installation is complete.

-------------------------
UNC Longleaf installation
-------------------------

Log into Longleaf using ssh and run the following code from the command line.

.. code-block::
   :caption: Downloading the github repository

   cd $HOME
   git clone https://github.com/Weeks-UNC/RNAvigate.git


Add the following line to .bash_profile or .bash_{ONYEN}.

.. code-block::
   :caption: Adding the package to PYTHONPATH

   export PYTHONPATH="$PYTHONPATH:$HOME/RNAvigate/"

Next we need to make sure that you have an Anaconda environment that includes
all dependencies and that that environment is available to Jupyter. This does
not change your default modules. They will be restored next time you log in.

.. code-block::
   :caption: Creating a Jupyter Notebook conda environment

   module rm python pymol pyrosetta
   module load anaconda/2019.10
   conda env create -f RNAvigate/env.yaml -n rnavigate
   source activate rnavigate
   python -m ipykernel install --user --name=rnavigate

If this occured without errors, exit longleaf.

OpenOnDemand is UNC's platform for using interactive programs within Longleaf.

1. Go to https://ondemand.rc.unc.edu/, log in, and start a jupyter notebook.

   1. Click on interactive apps, under servers, click on Jupyter Notebook.
   2. Enter the number of hours you will be working, 1 CPU, other fields blank.
   3. DONT FORGET TO SAVE BEFORE YOUR NOTEBOOK SHUTS DOWN.

2. Using the file navigation tree on the left, navigate to
   docs/source/examples/ and open rnavigate-example.ipynb.
3. In the upper right corner, click on "Python 3" and change it to "rnavigate"
   from the menu.
4. To test that this works, restart the kernel and run all cells.
5. On the JupyterLab home page, there should be an option to start a new python
   notebook using the "rnavigate" environment.
