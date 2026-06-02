FROM quay.io/jupyter/base-notebook:latest
ARG VERSION
RUN pip install --no-cache-dir "RNAvigate==${VERSION}"
