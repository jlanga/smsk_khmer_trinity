# smsk_khmer_trinity: a simple workflow for transcriptome assembly

[![Build Status](https://travis-ci.org/jlanga/smsk_khmer_trinity.svg?branch=master)](https://travis-ci.org/jlanga/smsk_khmer_trinity)

## 1. Description

This is a workflow for _de novo_ transcriptome assembly with Illumina reads. It

1. Trims reads with `Trimmomatic`

2. Performs digital normalization with `khmer`

3. Assembles with `trinity`

## 2. First steps

Just follow what is inside the `.travis.yml`

1. Update your system (Ubuntu Trusty in this case):
    ```sh
    sudo apt-get -qq update
    sudo apt-get install -y build-essential curl git
    ```

2. Clone this repo and get inside
    ```sh
    git clone https://github.com/jlanga/smsk_khmer_trinity
    cd smsk_khmer_trinity
    ```

3. "Activate" an environment (extend paths)
    ```sh
    source bin/activate
    ```

4. Install the required software

    ```sh
    bash scripts/install/brew.sh         # Brew itself
    bash scripts/install/from_brew.sh    # Custom pieces of software so you don't need to sudo
    bash scripts/install/from_pip3.sh    # Python3 packages
    bash scripts/install/from_tarball.sh # Software that cannot be installed with the previous methods
    ```

4. Execute the pipeline with test data:

    ```sh
    snakemake -j 16
    ```



## 3. File organization

The hierarchy of the folder is the one described in [A Quick Guide to Organizing Computational Biology Projects](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424):

```
smsk_khmer_trinity
├── .linuxbrew: brew files
├── bin: binaries and activate files.
├── data: raw data, hopefully links to backup data.
├── doc: logs, reports and figures.
├── README.md
├── results: processed data.
├── scripts: snakefiles, installing scriptos, python, R, etc scripts to process data.
└── src: additional source code, tarballs, etc.
```



## 4. Analyzing your data

"Just" edit the `config.yaml` with the paths to your fastq files and change parameters. In the section `diginorm_params` \ `max_table_size` type `4e9` because it's anoyingly slow to do tests with 16Gb of RAM.

Also raise Trinity's maximum memory usage if you need it.



## Links, References and Bibliography

- [A Quick Guide to Organizing Computational Biology Projects](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)

- [Snakemake—a scalable bioinformatics workflow engine](http://bioinformatics.oxfordjournals.org/content/28/19/2520)

- [Trimmomatic](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/)

- [khmer](https://khmer-protocols.readthedocs.io/en/latest/mrnaseq/)

- [trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
