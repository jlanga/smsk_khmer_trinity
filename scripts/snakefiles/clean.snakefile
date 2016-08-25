rule clean:
    shell:
        "rm -rf {RAW_DIR} {RAW_DOC}"
rule clean_raw:
    shell:
        "rm -rf {RAW_DIR} {RAW_DOC}"
