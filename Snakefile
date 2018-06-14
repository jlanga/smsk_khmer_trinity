shell.prefix("set -euo pipefail;")

configfile: "params.yml"
params = config.copy()

configfile: "samples.yml"
samples = config.copy()

del config

singularity: "docker://continuumio/miniconda3:4.4.10"

SAMPLES_PE = samples["samples_pe"] if samples["samples_pe"] is not None else []
SAMPLES_SE = samples["samples_se"] if samples["samples_se"] is not None else []

SAMPLES = [x for x in SAMPLES_PE] + [x for x in SAMPLES_SE]
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
