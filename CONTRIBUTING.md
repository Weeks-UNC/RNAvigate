# Contributing to RNAvigate

Thanks for taking the time to contribute! This guide covers how to set up a
development environment, how the test/lint/docs tooling fits together, and
what to expect when opening a pull request.

## Questions, bugs, and feature requests

Use [GitHub issues](https://github.com/Weeks-UNC/RNAvigate/issues) liberally
— for bug reports, unclear errors, feature requests (new file formats, new
visualizations), or just questions. You can also email the maintainer
directly (see the [README](https://github.com/Weeks-UNC/RNAvigate/blob/master/README.md)).

## Development environment setup

### Option 1: VS Code Dev Container (recommended)

This repo ships a [Dev Containers](https://containers.dev/) configuration
that pre-installs everything needed to write code, run tests, lint, and
build docs — no manual dependency management required.

Prerequisites:

- [Git](https://git-scm.com/downloads)
- [Visual Studio Code](https://code.visualstudio.com/)
- [Docker Desktop](https://www.docker.com/products/docker-desktop/)
- The **Dev Containers** VS Code extension (`ms-vscode-remote.remote-containers`)

Steps:

```bash
git clone https://github.com/Weeks-UNC/RNAvigate.git
code RNAvigate
```

VS Code will detect the `.devcontainer/` configuration and prompt **Reopen in
Container** in the bottom-right notification (or run **Dev Containers: Reopen
in Container** from the command palette). The first build takes a few
minutes; subsequent opens are fast since the image is cached.

The container is built on Python 3.12. All dependencies are pre-installed
into a virtual environment at `/opt/venv` (activated automatically in every
terminal), RNAvigate is installed in editable mode, and `pre-commit install`
runs automatically via `postCreateCommand` — so hooks are active immediately.

### Option 2: Manual setup

If you're not using the Dev Container, you'll need Python >= 3.9 and
[`uv`](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/Weeks-UNC/RNAvigate.git
cd RNAvigate
uv sync --extra dev --extra docs
uv run pre-commit install
```

`uv sync` creates a `.venv` and installs RNAvigate in editable mode along
with the `dev` extra (`pytest`, `pytest-cov`, `ruff`, `pre-commit`) and the
`docs` extra (`sphinx`, `sphinx-autobuild`, `sphinx-rtd-theme`, `nbsphinx`,
`myst-parser`, `ipykernel`). Prefix commands below with `uv run`, or
`source .venv/bin/activate` first.

Building the docs also requires [Pandoc](https://pandoc.org/installing.html)
(used by `nbsphinx` to render notebooks) — install it via your system package
manager (e.g. `apt-get install pandoc`, `brew install pandoc`).

> **Note:** the `pytest` and `sphinx` [pre-commit](https://github.com/Weeks-UNC/RNAvigate/blob/master/.pre-commit-config.yaml)
> hooks invoke `/opt/venv/bin/pytest` and `/opt/venv/bin/sphinx-build`
> directly, which only exists inside the Dev Container. In a manual setup
> those two hooks will fail to find their executable — run `pytest` and
> `sphinx-build` manually instead (see below) until this is made portable.

## Tools in the development environment

### Dependency management (`uv`)

Dependencies are declared in `pyproject.toml` and pinned in `uv.lock`. To add
or update a dependency, edit `pyproject.toml` and run `uv lock`, then
`uv sync --extra dev --extra docs` to apply it locally.

### Testing

```bash
pytest
```

Tests live under `tests/`, mirroring the package structure
(`tests/data/`, `tests/plots/`, `tests/analysis/`), and use fixtures from
`rnavigate.examples` for real data rather than synthetic fixtures. Coverage
reports:

```bash
pytest --cov=rnavigate --cov-report=html
```

### Linting and formatting

[Ruff](https://docs.astral.sh/ruff/) handles both linting and formatting:

```bash
ruff check --fix .
ruff format .
```

### Pre-commit hooks

Hooks run automatically on each `git commit` (see
[`.pre-commit-config.yaml`](https://github.com/Weeks-UNC/RNAvigate/blob/master/.pre-commit-config.yaml)):

- **ruff** / **ruff-format** — lint and format Python code
- **uv-lock** — keeps `uv.lock` in sync with `pyproject.toml`
- **nbstripout** — strips output cells from notebooks before committing
- **pytest** — runs the full test suite
- **sphinx** — builds the docs when files under `docs/` change

Run all hooks on demand without committing:

```bash
pre-commit run --all-files
```

### Building the docs

The documentation site is built with [Sphinx](https://www.sphinx-doc.org/).
Use `sphinx-autobuild` for a live-reloading local preview:

```bash
cd docs/
make livehtml
```

Once it reports `Serving on http://127.0.0.1:8000`, open that URL in a
browser. Changes to `.rst`, `.md`, or `.ipynb` files under `docs/source/`
trigger an automatic rebuild. Press `Ctrl+C` to stop the server.

For a one-off build (matching what CI runs):

```bash
sphinx-build docs/source docs/build/html
```

## Code style

See the full [style guide](https://rnavigate.readthedocs.io/en/latest/resources/style_guide.html)
for naming conventions, docstring format (NumPy-style), commit message
conventions, and line-length rules. The short version: `ruff` enforces
formatting and import order, PEP 8 naming, 88-character lines, and NumPy-style
docstrings on public functions/classes.

## Making changes / submitting a pull request

1. Create a branch off `master`.
2. Make the smallest change that accomplishes the goal — avoid unrelated
   refactoring in the same PR.
3. Add or update tests for any new or changed behavior.
4. If the change is structural (new module, moved class, new dependency),
   update [`ARCHITECTURE.md`](https://github.com/Weeks-UNC/RNAvigate/blob/master/ARCHITECTURE.md). If it's user-facing, add an
   entry to [`CHANGELOG.md`](https://github.com/Weeks-UNC/RNAvigate/blob/master/CHANGELOG.md)
   and check off the relevant item in
   [`TODO.md`](https://github.com/Weeks-UNC/RNAvigate/blob/master/TODO.md) if one exists.
5. Before opening the PR, make sure the full CI suite would pass locally:
   `ruff check .`, `ruff format --check .`, `pytest`, and (if docs changed)
   `sphinx-build docs/source docs/build/html`.
6. Open a pull request against `master` describing what changed and why.

For architectural context (module boundaries, data flow, class hierarchy)
before making non-trivial changes, see [`ARCHITECTURE.md`](https://github.com/Weeks-UNC/RNAvigate/blob/master/ARCHITECTURE.md).
If you're using an AI coding assistant, [`AGENTS.md`](https://github.com/Weeks-UNC/RNAvigate/blob/master/AGENTS.md) documents
the same conventions in a form meant for that workflow.

## Releasing a new version

Releases are currently triggered manually. Here are the steps:

1. Bump `version` in `pyproject.toml`.
2. In `CHANGELOG.md`, rename `## X.Y.Z (Coming soon)` to `## X.Y.Z (Month Day, Year)`.
3. Commit these changes (e.g. `chore: bump version to X.Y.Z`).
4. Tag the commit and push the tag:

   ```bash
   git tag vX.Y.Z
   git push origin vX.Y.Z
   ```

Pushing the tag triggers `.github/workflows/cd.yml`, which:

1. Reruns the CI suite (linting, testing, and documentation building).
2. Builds/publishes to PyPI using trusted publishing (no credentials required).
3. Waits for the new PyPI tag.
4. Builds/publishes to Docker Hub with 2 tags: `X.Y.Z` and `latest`.
5. Bioconda automatically detects the new tag and creates a PR.

If CI fails on the tagged commit, the PyPI publish and Docker build are
skipped — fix the issue, then delete and re-push the tag (or bump to the
next patch version and tag again).
