Remote HPC cluster via VS Code
==============================

This guide explains how to use a locally installed VS Code to create Jupyter notebooks
on a remote HPC cluster over SSH. You will need:

- An account on the cluster
- VS Code installed on your local machine
- An OpenSSH-compatible SSH client (built into macOS, Linux, and Windows 10+)

There are two types of connections you will manage:

- **Login node** — for editing files and lightweight tasks. Running
  computationally intensive work here is an abuse of shared resources.
- **Compute node** — for running actual calculations. Access requires
  submitting a job through the cluster's scheduler (typically SLURM).

.. contents:: Contents
   :local:


Install the Remote - SSH extension
------------------------------------

In VS Code, open the Extensions panel (**Ctrl+Shift+X** / **Cmd+Shift+X**),
search for ``Remote - SSH``, and install the extension by Microsoft
(``ms-vscode-remote.remote-ssh``). After installation a **Remote** icon
(two arrows) appears in the bottom-left corner of the VS Code window.


Set up SSH key authentication
------------------------------

Password authentication requires re-entering your credentials for every
VS Code operation. SSH key authentication avoids this.

1. Generate a key pair on your local machine:

.. code-block:: bash

   ssh-keygen -t rsa -b 2048

   Press **Enter** to accept the default file location (``~/.ssh/id_rsa``).
   A passphrase is optional.

2. Copy the public key to the cluster. On the cluster, append the contents
   of your local ``~/.ssh/id_rsa.pub`` to ``~/.ssh/authorized_keys``:

.. code-block:: bash

   # on your local machine
   cat ~/.ssh/id_rsa.pub

   # on the cluster (paste the output into authorized_keys)
   mkdir -p ~/.ssh
   touch ~/.ssh/authorized_keys
   # open authorized_keys in an editor and paste the public key on a new line

   Each key in ``authorized_keys`` must be on its own line with no line breaks
   within a single key entry.

.. warning::

   Never share your private key (``id_rsa``). Only the public key (``id_rsa.pub``)
   goes on the cluster.


Connect to the login node
--------------------------

1. In VS Code, click the **Remote** icon in the bottom-left corner and choose
   **Open SSH Configuration File…**, then select the config file under your
   home directory (``~/.ssh/config`` on macOS/Linux,
   ``C:\Users\<you>\.ssh\config`` on Windows).

2. Add an entry for the cluster login node:

.. code-block:: text

   Host <cluster>
       HostName <login-node>
       User <username>
       IdentityFile ~/.ssh/id_rsa

   Replace ``<cluster>`` with a short alias (e.g. ``mycluster``),
   ``<login-node>`` with the full hostname (e.g. ``login.cluster.edu``),
   and ``<username>`` with your cluster username. Save the file.

3. Click the **Remote** icon and choose **Connect to Host…**, then select
   ``<cluster>``. If prompted, confirm the platform is **Linux** and accept
   the host key fingerprint.

A new VS Code window opens connected to the login node. File editing and
terminal commands run on the cluster from here. Restrict this session to
lightweight tasks — editing, file management, and short tests only.

4. To disconnect: click the **Remote** icon in the bottom-left corner and
   choose **Close Remote Connection**.


Connect to a compute node
--------------------------

This section requires the login node connection from the previous section.

Step 1 — Update SSH config for compute nodes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compute nodes are typically not directly reachable from outside the cluster.
The SSH config below routes connections through the login node as a jump host.
The ``Host`` pattern must match your cluster's compute node naming convention
— check your cluster documentation for the correct pattern.

.. code-block:: text

   Host <cluster>
       HostName <login-node>

   Host *.<cluster-domain>
       ForwardAgent yes
       ForwardX11 yes
       ForwardX11Trusted yes
       IdentityFile ~/.ssh/id_rsa
       User <username>

   Host <compute-node-pattern>
       HostName %h
       ProxyJump <login-node>
       User <username>

``<compute-node-pattern>`` is a wildcard that matches your cluster's compute
node hostnames. For example, if compute nodes are named ``n001``, ``n002``,
etc., use ``n*``. If they are named ``gpu01.cluster.edu``, use
``*.cluster.edu``. Check your cluster's documentation for the correct pattern.

Step 2 — Request a compute node
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From the VS Code window connected to the login node, open a terminal
(**Terminal > New Terminal**) and submit an interactive SLURM job:

.. code-block:: bash

   srun --time=<hh:mm:ss> --partition=<interactive-partition> --nodes=1 --pty /bin/bash

Replace ``<hh:mm:ss>`` with the time you need and
``<interactive-partition>`` with your cluster's interactive partition name.
Once the job starts, a shell prompt on the compute node appears. Note the
compute node hostname:

.. code-block:: bash

   hostname
   # example output: n042.cluster.edu

You will need this hostname in the next step.

Step 3 — Connect VS Code to the compute node
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From the **same login-node VS Code window**, open the command palette
(**Ctrl+Shift+P** / **Cmd+Shift+P**), type **Remote-SSH: Connect to Host…**,
and enter the compute node hostname noted above (e.g. ``n042.cluster.edu``).

A new VS Code window opens connected to the compute node. Code executed in
this window runs on the allocated compute resources.

.. note::

   Closing the login node VS Code window also closes the compute node connection.
   Ending the SLURM interactive job (``exit`` in the terminal) does the same.

When finished, click the **Remote** icon and choose **Close Remote Connection**.

Connect to a compute node directly
------------------------------------

If you already have an interactive SLURM job running from a separate SSH
session (started outside VS Code), you can skip the login node window
entirely. Open VS Code, click the **Remote** icon, choose **Connect to
Host…**, and type the compute node hostname directly.

.. note::

   The interactive job must remain alive for the VS Code connection to persist.
   Exiting the ``srun`` session that launched the job will terminate the
   VS Code connection to the compute node.
