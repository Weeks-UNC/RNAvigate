FROM mcr.microsoft.com/vscode/devcontainers/python:3.9

COPY requirements.txt /tmp/pip-tmp/

RUN pip3 --disable-pip-version-check \
    --no-cache-dir install \
    -r /tmp/pip-tmp/requirements.txt \
    && rm -rf /tmp/pip-tmp

ENV PYTHONPATH "/workspaces/RNAvigate"

RUN apt-get update && apt-get install -y pandoc

# RUN mamba install --yes -c conda-forge -c bioconda --no-channel-priority \
#     'biopython' \
#     'openpyxl' \
#     'pandas' \
#     'py3dmol' \
#     'scikit-learn' \
#     'scipy' \
#     'seaborn' \
#     'mkdocs' \
#     'mkdocs-jupyter' \
#     'mkdocstrings' \
#     && mamba clean --all -f -y
