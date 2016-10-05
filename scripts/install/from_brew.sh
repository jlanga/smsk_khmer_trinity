brew update

# To make the pip3 and snakemake work
brew install \
    python3 \
    graphviz \

# Homebrew science packages to run the example
brew tap homebrew/science
brew update

# Bioinformatic tools
brew tap homebrew/science
brew install \
    homebrew/science/fastqc \
    homebrew/science/trimmomatic \
    homebrew/science/bowtie \


brew install \
    pigz
