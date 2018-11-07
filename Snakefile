import copy
import pandas as pd

shell.prefix("set -euo pipefail;")


configfile: "params.yml"
params = config.copy()
config = {}

configfile: "features.yml"
features = config.copy()
config = {}

del config

samples = pd.read_table("samples.tsv").set_index("sample")


singularity: "docker://continuumio/miniconda3:4.4.10"


SAMPLES_PE = samples[samples["type"] == "PE"].index.tolist()
SAMPLES_SE = samples[samples["type"] == "SE"].index.tolist()

SAMPLES = SAMPLES_PE + SAMPLES_SE
PAIRS = ["pepe", "pese"]

ALL_THREADS = 24


snakefiles = "src/snakefiles/"

include: snakefiles + "folders.py"
include: snakefiles + "clean.py"
include: snakefiles + "generic.py"
include: snakefiles + "raw.py"
include: snakefiles + "qc.py"
include: snakefiles + "diginorm.py"
include: snakefiles + "assembly.py"
include: snakefiles + "mapping.py"

rule all:
    """
    Run the entire pipeline:
        - Reports for raw, trimmed and normalized reads
        - Assembly
    """
    input:
        RAW + "multiqc_report.html",
        QC + "multiqc_report.html",
        NORM + "multiqc_report.html",
        ASSEMBLY + "Trinity.fasta",
        rules.mapping.input
