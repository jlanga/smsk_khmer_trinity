import os

rule fastqc:
    """FASTQC"""
    input:
        fq = "{filename}.fq.gz",
    output:
        html = protected("{filename}_fastqc.html"),
        zip = protected("{filename}_fastqc.zip")
    log:
        os.path.dirname("{prefix}.fq.gz") + "fastqc_{filename}.log"
    benchmark:
        os.path.dirname("{prefix}.fq.gz") + "fastqc_{filename}.bmk"
    conda:
        "generic.yml"
    shell:
        "fastqc --nogroup {input.fq} 2> {log} 1>&2"


rule fasta_fai:
    """Index a fasta"""
    input: "{whatever}.fasta"
    output: "{whatever}.fasta.fai"
    conda: "generic.yml"
    shell: "samtools faidx {input}"


rule multiqc_raw:
    '''MultiQC report over the FASTQC ones'''
    input:
        pe_files = expand(
            RAW + "{sample}_{pair}_fastqc.{extension}",
            sample=SAMPLES_PE,
            pair="1 2".split(),
            extension="html zip".split()
        ),
        se_files = expand(
            RAW + "{sample}_se_fastqc.{extension}",
            sample=SAMPLES_SE,
            extension="html zip".split()
        )
    output:
        html = protected(RAW + "multiqc_report.html")
    log:
        RAW + "multiqc.log"
    benchmark:
        RAW + "multiqc.bmk"
    conda:
        "generic.yml"
    shell:
        "multiqc --title Raw --filename {output.html} {RAW} 2> {log}"


rule multiqc_qc:
    input:
        pe_files = expand(
            QC + "{sample}_{pair}_fastqc.{extension}",
            sample=SAMPLES_PE,
            pair=PAIRS,
            extension="html zip".split()
        ),
        se_files = expand(
            QC + "{sample}_se_fastqc.{extension}",
            sample=SAMPLES_SE,
            extension="html zip".split()
        )
    output:
        html = protected(QC + "multiqc_report.html")
    threads:
        1
    params:
        folder = QC
    log:
        QC + "multiqc.log"
    benchmark:
        QC + "multiqc.bmk"
    conda:
        "generic.yml"
    shell:
        "multiqc --title QC --filename {output.html} {params.folder} 2> {log}"


rule multiquc_diginorm:
    input:
        files_pe = expand(
            NORM + "{sample}_{pair}_fastqc.{extension}",
            sample=SAMPLES_PE,
            pair=PAIRS,
            extension="html zip".split()
        ),
        files_se = expand(
            NORM + "{sample}_se_fastqc.{extension}",
            sample=SAMPLES_SE,
            extension="html zip".split()
        )
    output:
        html = protected(NORM + "multiqc_report.html")
    log:
        NORM + "multiqc.log"
    benchmark:
        NORM + "multiqc.bmk"
    conda:
        "generic.yml"
    shell:
        "multiqc --title Diginorm --filename {output.html} {NORM} 2> {log}"
