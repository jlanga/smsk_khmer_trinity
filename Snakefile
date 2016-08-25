shell.prefix("set -euo pipefail;")
configfile: "config.yaml"

SAMPLES_PE = config["samples_pe"]
SAMPLES_SE = config["samples_se"]

SAMPLES = [x for x in SAMPLES_PE] + [x for x in SAMPLES_SE]
PAIRS = ["pe_pe", "pe_se"]

BLOCK_THREADS = 99999

snakefiles = "scripts/snakefiles/"

include: snakefiles + "folders.snakefile"
include: snakefiles + "clean.snakefile"
include: snakefiles + "raw.snakefile"
include: snakefiles + "qc.snakefile"
include: snakefiles + "diginorm.snakefile"



rule all:
    input:
        expand(
            NORM_DIR + "{sample}.final.{pair}.fq.gz",
            sample = SAMPLES_PE,
            pair = PAIRS
        ),
        expand(
            NORM_DIR + "{sample}.final.se.fq.gz",
            sample = SAMPLES_SE
        ),
        expand(
            NORM_DOC + "{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_PE,
            pair = "pe_pe pe_se".split(),
            extension = "html zip".split()
        ),
        expand(
            NORM_DOC + "{sample}.final.se_fastqc.{extension}",
            sample = SAMPLES_SE,
            extension = "html zip".split()
        )
