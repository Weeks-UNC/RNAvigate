Installing RNAvigate
====================

Choose one of the installation options below:

- `Docker installation`_ - for local use with Windows/iOS/Linux
- `Anaconda installation`_ - for remote use on HPC clusters
- `Manual installation`_
- `UNC Longleaf installation`_

Docker installation
-------------------

To install RNAvigate locally, this is the easiest option.

Step 1 - Install Docker Desktop
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`Docker Desktop <https://www.docker.com/products/docker-desktop/>`_

Step 2 - Run RNAvigate docker image
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Open Docker Desktop. In the top search bar, type "psirving/rnavigate", and
click "run". Expand the "Optional settings" drop-down menu, and enter the
following options:

- Container name: RNAvigate
- Host port: 8888
- Volumes:

   - Host path: use the (...) button to choose a directory.
   - Container path: /home/jovyan/work

The choice of "Host path" is important. RNAvigate will *only* have access to this directory.
Choose a high level directory that includes data files you wish to work on, but not too high.
For example, your home directory is perhaps not appropriate.
A good choice for me would be my "weeks_lab" directory, where I do all of my lab related work.

Click "run" again. This will take you to a terminal, at the end of which you should see something like this::

   YYYY-MM-DD HH:MM:SS     To access the notebook, open this file in a browser:
   YYYY-MM-DD HH:MM:SS         file:///home/jovyan/.local/share/jupyter/runtime/nbserver-7-open.html
   YYYY-MM-DD HH:MM:SS     Or copy and paste one of these URLs:
   YYYY-MM-DD HH:MM:SS         http://12ab34de56fg:8888/?token=01234567890abcdefghijklmnopqrstuvwxyz01234567890
   YYYY-MM-DD HH:MM:SS      or http://127.0.0.1:8888/?token=01234567890abcdefghijklmnopqrstuvwxyz01234567890

Step 3 - Open Jupyter server
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Click on either the 2nd or 3rd link above, or copy/paste them into a browser.
This will open a Jupyter server.
Navigate to ``/home/jovyan/work/``, this is linked to the host path that you set up in step 2.
Your files are accessible here.
You can now open a new Jupyter Notebook and start exploring your data!

Step 4 - Save and quit
~~~~~~~~~~~~~~~~~~~~~~

Don't forget to save! To stop the session, go back to Docker Desktop, and click the stop button in the upper-right.

Step 5 - Restarting RNAvigate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To restart this session, go the "Containers" menu, click "start" on the RNAvigate container.
Next to that button, click (...) -> "view details" and go back to step 3.

Anaconda installation
---------------------

If you are working on shared HPC resources, the Docker method may not work. Use conda instead.

Step 1 - Download RNAvigate
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a terminal, navigate to the directory in which you'll save RNAvigate.
Clone the github repository.
Then, navigate into the new RNAvigate directory.

.. code-block:: bash

   git clone https://github.com/Weeks-UNC/RNAvigate.git RNAvigate_v1.0.0
   cd RNAvigate_v1.0.0

Step 2 - Create the conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a conda enviroment using the included environment.yml file.
This step can take a long time, be patient. Grab a coffee.

.. code-block:: bash

   conda env create -f environment.yml

Step 3 - Install RNAvigate as a Jupyter kernel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Activate the new conda environment.
Install the RNAvigate package in the new conda environment.
Then, install the new environment for Jupyter.

.. code-block:: bash

   conda activate RNAvigate_v1.0.0
   conda develop .
   python -m ipykernel install --user --name=RNAvigate_v1.0.0

Step 4 - Test the installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Test the installation by opening an example Jupyter Notebook.

You may need to follow additional guidance for using Jupyter Notebooks on HPCs.
This is dependent on the administration of your HPC, and outside the scope of this guide.

.. code-block:: bash

   jupyter notebook ./docs/examples/

Open a notebook and click "Kernel", "restart and run all".

If RNAvigate fails to import, add this code to the very top of the notebook.
Replace `/path/to/RNAvigate/` with the location of your RNAvigate directory.


.. code-block:: python

   import sys
   sys.path.append('/path/to/RNAvigate/')

This is an occassionally recurring issue with ``conda develop .`` used in step 3.
If you know of a solution, please let me know on Github issues.
I'm trying to avoid adding directly to system path,
so that different versions of RNAvigate work in different notebooks.

Manual installation
-------------------

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

Download the package from `Github <https://github.com/Weeks-UNC/RNAvigate>`_. Add
this directory your $PYTHONPATH, and you are good to go.

UNC Longleaf installation
-------------------------

Step 1 - Install RNAvigate
~~~~~~~~~~~~~~~~~~~~~~~~~~

First, we need to make sure that you have an Anaconda environment that includes
all dependencies and that that environment is available to Jupyter. This does
not change your default modules. They will be restored next time you log in.

.. code-block:: bash

   module rm python pymol pyrosetta
   module load anaconda/2019.10
   conda env create -f /proj/kweeks/bin/RNAvigate_v1.0.0/environment.yml
   source activate RNAvigate_v1.0.0
   conda develop /proj/kweeks/bin/RNAvigate_v1.0.0/
   python -m ipykernel install --user --name=RNAvigate_v1.0.0

If this occured without errors, exit longleaf and go to
UNC's `OpenOnDemand Service <https://ondemand.rc.unc.edu/>`_.
OpenOnDemand is UNC's platform for using interactive programs within Longleaf.

Step 2 - Open Jupyter on OnDemand
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Log in using your ONYEN, and start a jupyter notebook.

1. Click on interactive apps, under servers, click on Jupyter Notebook.
2. Enter the number of hours you will be working, 1 CPU, other fields blank.
3. Don't forget to save before you run out of time!

Step 3 - Test the installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Navigate to your data directory, and open a new notebook using the "RNAvigate_v1.0.0" option.
Type ``import rnavigate as rnav`` into the first code cell and run it.
If no errors occur, you're ready to starting exploring your data!
