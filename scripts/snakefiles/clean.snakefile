rule clean:
    shell:
        "rm -r "
            "{RAW_DIR} {RAW_DOC} "
            "{QC_DIR} {QC_DOC} "
            "{NORM_DIR} {NORM_DOC} "



rule clean_raw:
    shell:
        "rm -r {RAW_DIR} {RAW_DOC} "



rule clean_qc:
    shell:
        "rm -r {QC_DIR} {QC_DOC} "



rule clean_diginorm:
    shell:
        "rm -r {NORM_DIR} {NORM_DOC} "
