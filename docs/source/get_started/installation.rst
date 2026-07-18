Installing RNAvigate
====================

If you have any questions about installation, please submit a `GitHub Issue <https://github.com/Weeks-UNC/RNAvigate/issues>`_.

.. contents:: Choose the installation option that best fits your needs:
   :local:

- `Docker`_: beginner-friendly, local, any OS
- `Conda`_: command-line, local or remote, flexible but requires dependency management
- `Pip`_: command-line, local or remote, flexible but requires dependency management
- `VS Code`_: full IDE, local, worthwhile workflow improvement with more setup
- `UNC Longleaf`_: for UNC users working on the Longleaf HPC cluster
- :doc:`/resources/developer_installation`: VS Code devcontainer for contributing code
- :doc:`/resources/hpc_vscode`: local VS Code window with access to HPC files and compute resources

Each method follows the same three steps:

1. **Install dependencies**: set up RNAvigate and Jupyter in your environment
2. **Open a notebook**: launch JupyterLab in a browser or open a notebook in VS Code
3. **Test the installation**: run a sample plot to confirm everything works

=======================================================================================

Docker
---------------------------------------------------------------------------------------

With this method, you don't have to worry about managing dependencies or setting up
JupyterLab: the Docker container takes care of all that for you. The container runs a
JupyterLab server that you connect to from your browser.

Install Docker and Run the Container
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Install `Docker Desktop <https://www.docker.com/products/docker-desktop/>`_.

2. Open Docker Desktop and search for ``psirving/rnavigate``.

.. image:: ../images/docker_install/1_Docker_search.png

3. Make sure the tag is set to ``latest`` and click **Run**.

.. image:: ../images/docker_install/2_Docker_Run.png

4. Expand the **Optional settings** drop-down and enter the following:

   - **Container name**: ``RNAvigate``
   - **Host port**: ``8888``
   - **Volumes**:

     - **Host path**: use the (...) button to choose a directory.
     - **Container path**: ``/home/jovyan/work``

   .. note::

      Choose a **Host path** that contains your data files.
      RNAvigate will *only* have access to this directory.

.. image:: ../images/docker_install/3_Docker_setup.png

5. Click **Run**. A terminal window opens showing the container logs.

.. image:: ../images/docker_install/4_Docker_connect.png

Connect to JupyterLab from Docker
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

6. Click one of the links shown in the terminal to open JupyterLab in your browser.

7. Click on the **Work** directory; this is linked to the host path you set in step 4.

.. image:: ../images/docker_install/5_Jupyter_create.png

8. Open a new Jupyter Notebook.

Test Docker Installation
~~~~~~~~~~~~~~~~~~~~~~~~

See `Test Plot`_ below.

Stop and Restart
~~~~~~~~~~~~~~~~

- **Stop**: In Docker Desktop, click the **Stop** button on the container.
- **Restart**: Open Docker Desktop, go to **Containers**, click **Start** on the
  RNAvigate container, then click **(...) > View Details** and return to step 6.

=======================================================================================

Conda
---------------------------------------------------------------------------------------

Use this method if you prefer conda or mamba, or if Docker is unavailable. Care must be
taken to ensure environment dependencies are managed correctly, but this method is more
flexible than Docker and works on HPC clusters.

Install Conda and create the environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Download the ``environment.yml`` file.

.. code-block:: bash

   curl -O https://raw.githubusercontent.com/Weeks-UNC/RNAvigate/master/environment.yml

.. note::

   On Windows, open the URL above in a browser and save the file as ``environment.yml``.

2. Create the conda environment. This step can take several minutes: grab a coffee.

.. code-block:: bash

   conda env create -f environment.yml

.. note::

   To keep track of the RNAvigate version in the environment name, change the name in
   ``environment.yml`` before running the command above. For example, change
   ``name: rnavigate`` to ``name: rnavigate_v0.0.0``. Below, replace
   ``conda activate rnavigate`` with ``conda activate rnavigate_v0.0.0``.

3. Activate the environment.

.. code-block:: bash

   conda activate rnavigate

Connect to JupyterLab
~~~~~~~~~~~~~~~~~~~~~

4. Launch JupyterLab.

.. code-block:: bash

   jupyter lab

A browser window will open automatically. If not, copy the URL from the terminal.

Test the Installation
~~~~~~~~~~~~~~~~~~~~~

See `Test Plot`_ below.

=======================================================================================

Pip
---------------------------------------------------------------------------------------

Use this method if you prefer pip and virtual environments.

Install RNAvigate and Jupyter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Create and activate a virtual environment.

.. code-block:: bash

   python -m venv rnavigate-env
   source rnavigate-env/bin/activate

.. note::

   On Windows, use ``rnavigate-env\Scripts\activate`` instead.

2. Install RNAvigate and JupyterLab.

.. code-block:: bash

   pip install rnavigate jupyterlab

Connect to JupyterLab
~~~~~~~~~~~~~~~~~~~~~

3. Launch JupyterLab.

.. code-block:: bash

   jupyter lab

A browser window will open automatically. If not, copy/paste the URL from the terminal.

Test the Installation
~~~~~~~~~~~~~~~~~~~~~

See `Test Plot`_ below.

=======================================================================================

VS Code: more setup, local or remote, excellent interface
---------------------------------------------------------------------------------------

This method requires a little bit more setup than the others, but the experience is
worth it. VS Code runs Jupyter notebooks natively and can connect to generative AI for
code completion and debugging. If you are new to VS Code, this is a good way to get
familiar with it. If you are an experienced VS Code user, this workflow will be a big
improvement over the JupyterLab interface.

Install VS Code
~~~~~~~~~~~~~~~

1. Install `Visual Studio Code <https://code.visualstudio.com/>`_.

2. Open VS Code and install the following extensions from the Extensions panel
   (**Ctrl+Shift+X** / **Cmd+Shift+X**):

   - `Python <https://marketplace.visualstudio.com/items?itemName=ms-python.python>`_ (``ms-python.python``)
   - `Jupyter <https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter>`_ (``ms-toolsai.jupyter``)

3. Create an environment with pip or conda to install RNAvigate. 
   Either method works: VS Code detects pip and conda environments automatically.

.. code-block:: bash

   python -m venv rnavigate-env
   source rnavigate-env/bin/activate
   pip install rnavigate

.. note::

   On Windows, use ``rnavigate-env\Scripts\activate`` instead.

or

.. code-block:: bash

   conda env create -f environment.yml
   conda activate rnavigate

Open a Notebook
~~~~~~~~~~~~~~~

4. Use VS Code to open your working directory.

5. Create a new Jupyter Notebook: press **Ctrl+Shift+P** (**Cmd+Shift+P** on macOS),
   type **Jupyter: Create New Jupyter Notebook**, and press **Enter**.

6. Click **Select Kernel** in the top-right corner of the notebook, choose
   **Python Environments**, and select ``rnavigate-env``.

Test the Installation
~~~~~~~~~~~~~~~~~~~~~

See `Test Plot`_ below.

UNC Longleaf
------------

For UNC users working on the Longleaf HPC cluster.

Install the Environment
~~~~~~~~~~~~~~~~~~~~~~~

Load the required modules and create a conda environment.
This does not change your default modules: they will be restored on your next login.

.. code-block:: bash

   module rm python pymol pyrosetta
   module load anaconda/2019.10

Follow the instructions in the `Install Conda and create the environment`_ section
above (don't run ``jupyter lab``) to create an environment with RNAvigate.

If this completes without errors, exit Longleaf and open UNC's
`OpenOnDemand Service <https://ondemand.rc.unc.edu/>`_.

Connect to Jupyter on OnDemand
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Log in with your ONYEN and start a Jupyter Notebook.

1. Click **Interactive Apps** and under **Servers** click **Jupyter Notebook**.
2. Enter the number of hours you will need, set CPU to 1, and leave other fields blank.
3. Don't forget to save your work before time runs out!

Test the Installation
~~~~~~~~~~~~~~~~~~~~~

See `Test Plot`_ below.

=======================================================================================

Test Plot
---------

Run the following in a Jupyter Notebook cell to confirm RNAvigate is working correctly.
This loads a built-in example dataset and generates an arc plot of DMS-MaP reactivity.

.. code-block:: python

   import rnavigate as rnav
   from rnavigate.examples import tpp

   rnav.plot_arcs(
       samples=[tpp],
       sequence="ss",
       structure="ss",
       profile="dmsmap",
   )

If the plot appears without errors, your installation is working correctly.
