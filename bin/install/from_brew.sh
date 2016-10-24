brew update

# To make the pip3 and snakemake work
brew install \
    python3 \
    graphviz \

# To make trinity work
brew install \
	python


# Homebrew science packages to run the example
brew tap homebrew/science
brew update

# Bioinformatic tools
brew tap homebrew/science
brew install \
    homebrew/science/fastqc \
    homebrew/science/trimmomatic \
    homebrew/science/bowtie \
    homebrew/science/trinity --without-express 


brew install \
    pigz
