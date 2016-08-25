shell.prefix("set -euo pipefail;")
configfile: "config.yaml"

SAMPLES_PE = config["samples_pe"]
SAMPLES_SE = config["samples_se"]

SAMPLES = [x for x in SAMPLES_PE] + [x for x in SAMPLES_SE]

snakefiles = "scripts/snakefiles/"

include: snakefiles + "folders.snakefile"
include: snakefiles + "clean.snakefile"
include: snakefiles + "raw.snakefile"

rule all:
    input:
        expand(
            RAW_DOC + "{sample}_{end}_fastqc.{extension}",
            sample = SAMPLES_PE,
            end = "1 2".split(),
            extension = ["html", "zip"]
        ),
        expand(
            RAW_DOC + "{sample}_se_fastqc.{extension}",
            sample = SAMPLES_SE,
            extension = ["html", "zip"]
        )
