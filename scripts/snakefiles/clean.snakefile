rule clean:
    shell:
        "rm -rf {RAW_DIR} {RAW_DOC} {QC_DIR} {QC_DOC}"

rule clean_raw:
    shell:
        "rm -rf {RAW_DIR} {RAW_DOC}"

rule clean_qc:
    shell:
        "rm -rf {QC_DIR} {QC_DOC}"
