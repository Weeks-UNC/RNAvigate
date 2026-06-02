FROM mcr.microsoft.com/devcontainers/python:3.12

# Install UV from the official image
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

# System deps (pandoc required by nbsphinx for docs builds)
RUN apt-get update && apt-get install -y pandoc && rm -rf /var/lib/apt/lists/*

# Pin the venv to a stable path outside the workspace so it survives rebuilds
# and is not affected by the workspace mount overwriting it.
ENV UV_PROJECT_ENVIRONMENT=/opt/venv
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Pre-install all Python dependencies at image-build time using the lockfile.
# --no-install-project: skips installing the project itself because the source
# is not present during the build — it gets mounted at runtime.
# The project is installed in editable mode by postCreateCommand (see devcontainer.json).
COPY pyproject.toml uv.lock /tmp/rnavigate-build/
RUN cd /tmp/rnavigate-build \
    && mkdir -p rnavigate && touch rnavigate/__init__.py \
    && uv sync --frozen --extra dev --extra docs --no-install-project \
    && rm -rf /tmp/rnavigate-build \
    && chown -R vscode:vscode /opt/venv
