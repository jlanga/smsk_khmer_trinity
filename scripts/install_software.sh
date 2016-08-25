#!/usr/bin/env bash
set -euo pipefail

## 1. Pip packages
pip3 install \
    snakemake \
    pyyaml \
    docutils



## 2. Homebrew packages
brew update

brew tap homebrew/science
brew install \
    homebrew/science/fastqc \
    homebrew/science/trimmomatic \
    homebrew/science/bowtie \
    pigz




## 3. Custom packages
## wget tar.gz; tar xvf tar.gz; make; make test; make install


mkdir -p src/
pushd src/

### Screed
echo "Installing screed"
git clone https://github.com/ged-lab/screed.git
pushd screed/
git checkout protocols-v0.8.3
python setup.py install
popd



### Khmer
echo "Installing khmer"
git clone https://github.com/ged-lab/khmer.git
pushd khmer/
git checkout protocols-v0.8.3
make
make install
popd



### Installing eel-pond scripts
# eel-pond scripts
git clone https://github.com/ctb/eel-pond.git
pushd eel-pond/
git checkout protocols-v0.8.3
popd



# Trinity-2.1.1 - Impossible to install via brew. issues with express and
# protobuf: https://github.com/Homebrew/homebrew-science/issues/3828
echo "Installing Trinity-2.2.0"
curl \
    --location \
    --output trinityrnaseq-2.2.0.tar.gz \
    https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.2.0.tar.gz
tar xzvf trinityrnaseq-2.2.0.tar.gz
pushd trinityrnaseq-2.2.0/
make -j
make -j plugins
pushd sample_data/test_Trinity_Assembly/
make -j
popd
ln -s ../src/trinityrnaseq-2.2.0/Trinity ../../bin/Trinity
popd


popd
