#!/usr/bin/env bash

set -euxo pipefail

# We do this conditionally because it saves us some downloading if the
# version is the same.
if [ ! -d "$HOME"/miniconda/bin ]; then
    echo "Downloading and installing conda"
    wget \
        --continue \
        https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        --output-document=miniconda.sh
    rm -rf "$HOME"/miniconda
    bash miniconda.sh -b -p "$HOME"/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    hash -r
    conda config --set always_yes yes --set changeps1 no
else
    echo "conda already installed. updating conda"
    export PATH="$HOME/miniconda/bin:$PATH"
    conda update -q --yes conda
fi

# Useful for debugging any issues with conda
conda info -a

# Use bioconda and conda-forge
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install --yes snakemake=5.0.0 yamllint=1.11.1 pylint=1.9.2
conda clean --all --yes
pip install pathspec
