Developer installation
======================

This guide sets up a full development environment using VS Code's Dev Containers
feature. Everything needed to write code, run tests, lint, and build docs is
pre-installed in the container: no manual dependency management required.

.. contents:: Contents
   :local:

Prerequisites
-------------

Install the following before getting started:

- `Git <https://git-scm.com/downloads>`_
- `Visual Studio Code <https://code.visualstudio.com/>`_
- `Docker Desktop <https://www.docker.com/products/docker-desktop/>`_
- The **Dev Containers** extension for VS Code (``ms-vscode-remote.remote-containers``)

Clone the repository
--------------------

.. code-block:: bash

   git clone https://github.com/Weeks-UNC/RNAvigate.git


Open in VS Code
---------------

Open the cloned directory in VS Code:

.. code-block:: bash

   code RNAvigate

Or from within VS Code: **File > Open Folder** and select the ``RNAvigate`` directory.


Reopen in Dev Container
-----------------------

VS Code will detect the ``.devcontainer`` configuration and prompt you
to reopen in the container.

1. Click **Reopen in Container** in the notification that appears in the
   bottom-right corner of VS Code.

   Alternatively, press **Ctrl+Shift+P** (**Cmd+Shift+P** on macOS),
   type **Dev Containers: Reopen in Container**, and press **Enter**.

2. The first time you open the container, Docker will build the image and
   install all dependencies. This may take a few minutes. Subsequent opens
   are much faster because the image is cached.

3. Once the container is ready, VS Code reloads with the full development
   environment active and all extensions installed.


Tools in the Dev Container
--------------------------

The dev container is built on Python 3.12 and includes the following tools.

Python Environment
~~~~~~~~~~~~~~~~~~

All dependencies are pre-installed in a virtual environment at ``/opt/venv``
and activated automatically in every terminal. RNAvigate itself is installed
in editable mode, so changes to the source code take effect immediately
without reinstalling.

uv
~~

`uv <https://docs.astral.sh/uv/>`_ is the package manager for the project.
Dependencies are declared in ``pyproject.toml`` and pinned in ``uv.lock``.
To add or update a dependency, edit ``pyproject.toml`` and run:

.. code-block:: bash

   uv sync --extra dev --extra docs

The ``dev`` extra includes ``pytest``, ``pytest-cov``, ``ruff``, and ``pre-commit``.
The ``docs`` extra includes ``sphinx``, ``sphinx-autobuild``, ``sphinx-rtd-theme``,
``nbsphinx``, and ``ipykernel``.

Testing
~~~~~~~

Run the test suite with pytest:

.. code-block:: bash

   pytest

Coverage reports can be generated with:

.. code-block:: bash

   pytest --cov=rnavigate --cov-report=html

Linting and Formatting
~~~~~~~~~~~~~~~~~~~~~~

`Ruff <https://docs.astral.sh/ruff/>`_ handles both linting and formatting,
and is configured to format Python files on save in VS Code.
To run manually:

.. code-block:: bash

   ruff check --fix .
   ruff format .

Pre-commit Hooks
~~~~~~~~~~~~~~~~

Pre-commit hooks run automatically on each ``git commit`` and enforce the
following checks:

- **ruff**: lints Python code and applies safe auto-fixes
- **ruff-format**: formats Python code
- **uv-lock**: keeps ``uv.lock`` in sync with ``pyproject.toml``
- **nbstripout**: strips output cells from notebooks before committing
- **pytest**: runs the full test suite

To run all hooks manually without committing:

.. code-block:: bash

   pre-commit run --all-files

Documentation
~~~~~~~~~~~~~

This documentation website is built with `Sphinx <https://www.sphinx-doc.org/>`_.
You can use ``sphinx-autobuild`` from the VS Code terminal to preview it locally.

.. code-block:: bash

   cd docs/
   make livehtml

If the build is successful, you will see a message like:

.. code-block:: text

   [sphinx-autobuild] Building...
   ...
   ...lots of messages about the build process...
   ...
   build succeeded.
   The HTML pages are in build.
   [sphinx-autobuild] Serving on http://127.0.0.1:8000

``Alt + click`` on the URL to open the docs in your browser.
Changes to ``.rst`` or ``.ipynb`` files trigger an automatic rebuild.
Type ``Ctrl + C`` in the terminal to stop the server when you're done.

VS Code Extensions
~~~~~~~~~~~~~~~~~~

The following extensions are installed automatically in the container:

- **Python** (``ms-python.python``): Python language support
- **Pylance** (``ms-python.vscode-pylance``): type checking and IntelliSense
- **Ruff** (``charliermarsh.ruff``): inline linting and format-on-save
- **Jupyter** (``ms-toolsai.jupyter``): Jupyter notebook support
- **Better Comments** (``aaron-bond.better-comments``): styled code annotations
- **Markdown Preview Enhanced** (``shd101wyy.markdown-preview-enhanced``): live Markdown preview
- **Text Tables** (``RomanPeshkov.vscode-text-tables``): RST and Markdown table editing
- **VS Code Icons** (``vscode-icons-team.vscode-icons``): file icons
- **Claude Code** (``Anthropic.claude-code``): AI coding assistant
