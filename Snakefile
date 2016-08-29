shell.prefix("set -euo pipefail;")
configfile: "config.yaml"

SAMPLES_PE = config["samples_pe"]
SAMPLES_SE = config["samples_se"]

SAMPLES = [x for x in SAMPLES_PE] + [x for x in SAMPLES_SE]
PAIRS = ["pe_pe", "pe_se"]

BLOCK_THREADS = 99999
ALL_THREADS = 24

snakefiles = "scripts/snakefiles/"

include: snakefiles + "folders.snakefile"
include: snakefiles + "clean.snakefile"
include: snakefiles + "raw.snakefile"
include: snakefiles + "qc.snakefile"
include: snakefiles + "diginorm.snakefile"
include: snakefiles + "assembly.snakefile"


rule all:
    input:
        ASSEMBLY_DIR + "Trinity.fasta",
        RAW_DOC + "multiqc_report.html",
        QC_DOC + "multiqc_report.html",
        NORM_DOC + "multiqc_report.html"
