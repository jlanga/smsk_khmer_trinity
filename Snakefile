import copy

shell.prefix("set -euo pipefail;")


configfile: "params.yml"
params = config.copy()
config = {}

configfile: "samples.yml"
samples = config.copy()
config = {}

configfile: "features.yml"
features = config.copy()
config = {}

del config


singularity: "docker://continuumio/miniconda3:4.4.10"


SAMPLES_PE = [
    sample for sample in samples
    if samples[sample]["type"] == "PE"
]

SAMPLES_SE = [
    sample for sample in samples
    if samples[sample]["type"] == "SE"
]

SAMPLES = samples.keys()
PAIRS = ["pepe", "pese"]

ALL_THREADS = 64


snakefiles = "src/snakefiles/"

include: snakefiles + "folders.py"
include: snakefiles + "clean.py"
include: snakefiles + "generic.py"
include: snakefiles + "raw.py"
include: snakefiles + "qc.py"
include: snakefiles + "diginorm.py"
include: snakefiles + "assembly.py"

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
        ASSEMBLY + "Trinity.fasta"
