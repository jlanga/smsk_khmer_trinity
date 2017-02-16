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
include: snakefiles + "filtering"
include: snakefiles + "transrate"
include: snakefiles + "tissue"

rule all:
    input:
        transrate + "assemblies.csv",
        raw + "multiqc_report.html",
        qc + "multiqc_report.html",
        norm + "multiqc_report.html",
        expand(
            tissue + "ids_{sample}.tsv",
            sample = SAMPLES_PE
        ) 
