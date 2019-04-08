# smsk_khmer_trinity: a simple workflow for transcriptome assembly

[![Build Status](https://travis-ci.org/jlanga/smsk_khmer_trinity.svg?branch=master)](https://travis-ci.org/jlanga/smsk_khmer_trinity)

## 1. Description

This is a workflow for _de novo_ transcriptome assembly with Illumina reads. It

1. Trims reads with `Trimmomatic`

2. Performs digital normalization with `khmer`

3. Assembles with `trinity`

## 2. First steps

Just follow what is inside the `.travis.yml`

1. Install `conda`

2. Clone this repo

3. Add your samples to `config.yaml`

4. Run snakemake: `snakemake --use-conda -j`

## 3. File organization

The hierarchy of the folder is the one described in [A Quick Guide to Organizing Computational Biology Projects](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424):

```none
smsk_khmer_trinity
├── bin: binaries, scripts and environment files.
├── data: raw data, hopefully links to backup data.
├── README.md - This
├── results: processed data.
|    ├── raw: links to raw data
|    ├── qc: processed reads with trimmomatic
|    ├── diginorm: digital normalization
|    ├── assembly: Trinity output
|    ├── filtering: TPM per loci filtering
|    ├── tissue: per sample quantification
|    └── transrate: assembly and filtering statistics
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

- [kallisto](https://pachterlab.github.io/kallisto/)

- [sleuth](http://pachterlab.github.io/sleuth/)

- [transrate](hibberdlab.com/transrate/)
