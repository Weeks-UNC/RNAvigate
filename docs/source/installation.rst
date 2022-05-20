============
Installation
============

.. note::

   If you are using UNC's Longleaf computing cluster, skip ahead to "Longleaf
   installation".

Open a terminal and navigate to the directory you'll be saving StarMapper, then
clone the github repository. Navigate into the StarMapper directory and install
the package.

.. code-block::
   :caption: Downloading the github repository

   git clone https://github.com/Weeks-UNC/StarMapper.git
   cd StarMapper
   pip install .

If you are using Anaconda to manage environments, the easiest way to ensure an
appropriate environment is to create one using the included env.yaml file. If
not, open the yaml file in a text editor. There you will find the dependencies
that you will need to install.

.. code-block::
   :caption: Creating a Jupyter Notebook conda environment

   conda env create -f StarMapper/env.yaml -n starmapper
   source activate starmapper
   python -m ipykernel install --user --name=starmapper

If this occured without errors, your installation is complete.

-------------------------
UNC Longleaf installation
-------------------------

Log into Longleaf using ssh and run the following code from the command line.

.. code-block::
   :caption: Downloading the github repository

   cd $HOME
   git clone https://github.com/Weeks-UNC/StarMapper.git


Add the following line to .bash_profile or .bash_{ONYEN}.

.. code-block::
   :caption: Adding the package to PYTHONPATH

   export PYTHONPATH="$PYTHONPATH:$HOME/StarMapper/"

Next we need to make sure that you have an Anaconda environment that includes
all dependencies and that that environment is available to Jupyter. This does
not change your default modules. They will be restored next time you log in.

.. code-block::
   :caption: Creating a Jupyter Notebook conda environment

   module rm python pymol pyrosetta
   module load anaconda/2019.10
   conda env create -f StarMapper/env.yaml -n starmapper
   source activate starmapper
   python -m ipykernel install --user --name=starmapper

If this occured without errors, exit longleaf.

OpenOnDemand is UNC's platform for using interactive programs within Longleaf.

1. Go to https://ondemand.rc.unc.edu/, log in, and start a jupyter notebook.

   1. Click on interactive apps, under servers, click on Jupyter Notebook.
   2. Enter the number of hours you will be working, 1 CPU, other fields blank.
   3. DONT FORGET TO SAVE BEFORE YOUR NOTEBOOK SHUTS DOWN.

2. Using the file navigation tree on the left, navigate to
   docs/source/examples/ and open starmapper-example.ipynb.
3. In the upper right corner, click on "Python 3" and change it to "starmapper"
   from the menu.
4. To test that this works, restart the kernel and run all cells.
5. On the JupyterLab home page, there should be an option to start a new python
   notebook using the "starmapper" environment.
