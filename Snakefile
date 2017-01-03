shell.prefix("set -euo pipefail;")
configfile: "config.yaml"


SAMPLES_PE = config["samples_pe"] if config["samples_pe"] is not None else []
SAMPLES_SE = config["samples_se"] if config["samples_se"] is not None else []

SAMPLES = [x for x in SAMPLES_PE] + [x for x in SAMPLES_SE]
PAIRS = ["pe_pe", "pe_se"]

ALL_THREADS = 64

snakefiles = "bin/snakefiles/"

include: snakefiles + "folders"
include: snakefiles + "clean"
include: snakefiles + "raw"
include: snakefiles + "qc"
include: snakefiles + "diginorm"
include: snakefiles + "assembly"


rule all:
    input:
        ASSEMBLY_DIR + "Trinity.fasta",
        RAW_DIR + "multiqc_report.html",
        QC_DIR + "multiqc_report.html",
        NORM_DIR + "multiqc_report.html"
