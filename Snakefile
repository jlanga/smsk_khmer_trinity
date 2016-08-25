shell.prefix("set -euo pipefail;")
configfile: "config.yaml"


snakefiles = "scripts/snakefiles/"

include: snakefiles + "folders.snakefile"
include: snakefiles + "clean.snakefile"
include: snakefiles + "raw.snakefile"

rule all:
    input:
        expand(
            RAW_DOC + "{sample}_1_fastqc.html",
            sample = config["samples_pe"]
        ),
        expand(
            RAW_DOC + "{sample}_se_fastqc.html",
            sample = config["samples_se"]
        )
