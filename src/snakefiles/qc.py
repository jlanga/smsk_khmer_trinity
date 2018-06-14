def get_adaptor_pe(wildcards):
    """Get the adaptor from a sample"""
    return samples["samples_pe"][wildcards.sample]["adaptor"]


def get_phred_pe(wildcards):
    """Get the adaptor from a sample"""
    return samples["samples_pe"][wildcards.sample]["phred"]


rule qc_trimmomatic_pe:
    """
    Run trimmomatic on paired end mode to eliminate Illumina adaptors and
    remove low quality regions and reads.
    Inputs _1 and _2 are piped through gzip/pigz.
    Outputs _1 and _2 are piped to gzip/pigz (level 1, will be deleted).
    Outputs _3 and _4 are compressed with the builtin compressor from
    Trimmomatic. Further on they are catted and compressed with gzip/pigz
    (level 9).
    Note: The cut -f 1 -d " " is to remove additional fields in the FASTQ
    header. It is done posterior to the trimming since the output comes
    slower than the input is read.
    Number of threads used:
        4 for trimmomatic
        2 for gzip inputs
        2 for gzip outputs
        Total: 8
    """
    input:
        forward = RAW + "{sample}_1.fq.gz",
        reverse = RAW + "{sample}_2.fq.gz"
    output:
        forward = temp(QC + "{sample}_1.fq.gz"),
        reverse = temp(QC + "{sample}_2.fq.gz"),
        unpaired = protected(QC + "{sample}_pese.fq.gz")
    params:
        unpaired_1 = QC + "{sample}_3.fq.gz",
        unpaired_2 = QC + "{sample}_4.fq.gz",
        adaptor = get_adaptor_pe,
        phred = get_phred_pe,
        trimmomatic_params = params["trimmomatic"]["extra"]
    log:
        QC + "trimmomatic_pe_{sample}.log"
    benchmark:
        QC + "trimmomatic_pe_{sample}.bmk"
    priority:
        50
    threads:
        4
    conda:
        "qc.yml"
    shell:
        """
        trimmomatic PE \
            -threads {threads} \
            -{params.phred} \
            <(gzip --decompress --stdout {input.forward}) \
            <(gzip --decompress --stdout {input.reverse}) \
            >(  cut --fields 1 --delimiter \" \" \
                | pigz --fast > {output.forward}) \
            {params.unpaired_1} \
            >(  cut --fields 1 --delimiter \" \" \
                | pigz --fast > {output.reverse}) \
            {params.unpaired_2} \
            ILLUMINACLIP:{params.adaptor}:2:30:10 \
            {params.trimmomatic_params} 2> {log}

        gzip --decompress --stdout {params.unpaired_1} {params.unpaired_2} \
        | cut --fields 1 --delimiter \" \" \
        | pigz --best \
        > {output.unpaired} 2>> {log}

        rm {params.unpaired_1} {params.unpaired_2} 2>> {log}
        """


def get_adaptor_se(wildcards):
    """Get the adaptor from a sample"""
    return samples["samples_se"][wildcards.sample]["adaptor"]


def get_phred_se(wildcards):
    """Get the adaptor from a sample"""
    return samples["samples_se"][wildcards.sample]["phred"]


rule qc_trimmomatic_se:
    """
    Run trimmomatic on single end mode to eliminate Illumina adaptors and
        remove low quality regions and reads.
    Input is piped through gzip/pigz.
    Output is piped to gzip.
    Threads used:
        4 for trimmomatic
        1 for gzip input
        1 for gzip output
    """
    input:
        single = RAW + "{sample}_se.fq.gz",
    output:
        single = protected(QC + "{sample}_se.fq.gz")
    params:
        adaptor = get_adaptor_se,
        phred = get_phred_se,
        trimmomatic_params = params["trimmomatic"]["extra"]
    priority:
        50
    log:
        QC + "trimmomatic_se_{sample}.log"
    benchmark:
        QC + "trimmomatic_se_{sample}.bmk"
    threads:
        ALL_THREADS
    conda:
        "qc.yml"
    shell:
        """
        trimmomatic SE \
            -threads {threads} \
            -{params.phred} \
            <(gzip --decompress --stdout {input.single}) \
            >(  cut --fields 1 --delimiter \" \" \
                | pigz --best > {output.single}) \
            ILLUMINACLIP:{params.adaptor}:2:30:10 \
            {params.trimmomatic_params} \
        2> {log}
        """


rule qc_interleave_pepe:
    """
    From the adaptor free _1 and _2 , interleave the reads.
    Read the inputs, interleave, filter the stream and compress.
    """
    input:
        forward = QC + "{sample}_1.fq.gz",
        reverse = QC + "{sample}_2.fq.gz"
    output:
        interleaved = protected(QC + "{sample}_pepe.fq.gz")
    log:
        QC + "interleave_pe_{sample}.log"
    benchmark:
        QC + "interleave_pe_{sample}.bmk"
    conda:
        "qc.yml"
    shell:
        """
        (interleave-reads.py \
            <(gzip -dc {input.forward}) \
            <(gzip -dc {input.reverse}) \
        | pigz --best \
        > {output.interleaved}) \
        2> {log}
        """


rule qc_results:
    """Generate only the resulting files, not the reports"""
    input:
        pe_files = expand(
            QC + "{sample}_{pair}.fq.gz",
            sample=SAMPLES_PE,
            pair=PAIRS
        ),
        se_files = expand(
            QC + "{sample}_se.fq.gz",
            sample=SAMPLES_SE
        )


rule qc_doc:
    input:
        html = QC + "multiqc_report.html"


rule qc:
    """qc_results + qc_doc"""
    input:
        rules.qc_results.input,
        rules.qc_doc.input
