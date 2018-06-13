shell.prefix("set -euo pipefail;")

configfile: "config.yml"

singularity: "docker://continuumio/miniconda3:4.4.10"

SAMPLES_PE = config["samples_pe"] if config["samples_pe"] is not None else []
SAMPLES_SE = config["samples_se"] if config["samples_se"] is not None else []

SAMPLES = [x for x in SAMPLES_PE] + [x for x in SAMPLES_SE]
PAIRS = ["pe_pe", "pe_se"]

ALL_THREADS = 64

snakefiles = "src/snakefiles/"

include: snakefiles + "folders.py"
include: snakefiles + "clean.py"
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
        raw + "multiqc_report.html",
        qc + "multiqc_report.html",
        norm + "multiqc_report.html",
        assembly + "Trinity.fasta"
