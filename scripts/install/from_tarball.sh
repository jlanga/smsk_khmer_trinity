#!/usr/bin/env bash

set -euo pipefail # Bash unoficial strict mode


mkdir -p src/
pushd src/

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
