def get_adaptor_pe(wildcards):
    """Get the adaptor from a sample"""
    return config["samples_pe"][wildcards.sample]["adaptor"]


def get_phred_pe(wildcards):
    """Get the adaptor from a sample"""
    return config["samples_pe"][wildcards.sample]["phred"]


rule qc_trimmomatic_pe:
    """
    Run trimmomatic on paired end mode to eliminate Illumina adaptors and
    remove low quality regions and reads.
    Inputs _1 and _2 are piped through gzip/pigz.
    Outputs _1 and _2 are piped to gzip/pigz (level 9).
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
        unpaired = protected(QC + "{sample}.final.pe_se.fq.gz")
    params:
        unpaired_1 = QC + "{sample}_3.fq.gz",
        unpaired_2 = QC + "{sample}_4.fq.gz",
        adaptor = get_adaptor_pe,
        phred = get_phred_pe,
        trimmomatic_params = config["trimmomatic_params"]
    log:
        QC + "trimmomatic_pe_{sample}.log"
    benchmark:
        QC + "trimmomatic_pe_{sample}.bmk"
    priority:
        50
    threads:
        ALL_THREADS
    conda:
        "qc.yml"
    shell:
        """
        (trimmomatic PE \
            -threads {threads} \
            -{params.phred} \
            <(gzip --decompress --stdout {input.forward}) \
            <(gzip --decompress --stdout {input.reverse}) \
            >(  cut --fields 1 --delimiter \" \" \
                | pigz --best > {output.forward}) \
            {params.unpaired_1} \
            >(  cut --fields 1 --delimiter \" \" \
                | pigz --best > {output.reverse}) \
            {params.unpaired_2} \
            ILLUMINACLIP:{params.adaptor}:2:30:10 \
            {params.trimmomatic_params}
        gzip --decompress --stdout {params.unpaired_1} {params.unpaired_2} \
        | cut --fields 1 --delimiter \" \" \
        | pigz --best \
        > {output.unpaired}
        rm {params.unpaired_1} {params.unpaired_2}) \
        2> {log}
        """


def get_adaptor_se(wildcards):
    """Get the adaptor from a sample"""
    return config["samples_se"][wildcards.sample]["adaptor"]


def get_phred_se(wildcards):
    """Get the adaptor from a sample"""
    return config["samples_se"][wildcards.sample]["phred"]


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
        single = protected(QC + "{sample}.final.se.fq.gz")
    params:
        adaptor = get_adaptor_se,
        phred = get_phred_se,
        trimmomatic_params = config["trimmomatic_params"]
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


rule qc_decompress_pe_sample:
    input:
        forward = QC + "{sample}_1.fq.gz",
        reverse = QC + "{sample}_2.fq.gz"
    output:
        forward = temp(QC + "{sample}_1.fq"),
        reverse = temp(QC + "{sample}_2.fq")
    log:
        QC + "decompress_pe_{sample}.log"
    benchmark:
        QC + "decompress_pe_{sample}.bmk"
    conda:
        "qc.yml"
    shell:
        "gzip --decompress --keep {input.forward} {input.reverse} 2> {log}"


rule qc_interleave_pe_pe:
    """
    From the adaptor free _1 and _2 , interleave the reads.
    Read the inputs, interleave, filter the stream and compress.
    """
    input:
        forward = QC + "{sample}_1.fq.gz",
        reverse = QC + "{sample}_2.fq.gz"
    output:
        interleaved = protected(QC + "{sample}.final.pe_pe.fq.gz")
    log:
        QC + "interleave_pe_{sample}.log"
    benchmark:
        QC + "interleave_pe_{sample}.bmk"
    conda:
        "qc.yml"
    shell:
        """
        (interleave-reads.py \
            <(gzip --decompress --stdout {input.forward}) \
            <(gzip --decompress --stdout {input.reverse}) \
        | pigz --best \
        > {output.interleaved}) \
        2> {log}
        """


rule qc_results:
    """Generate only the resulting files, not the reports"""
    input:
        pe_files = expand(
            QC + "{sample}.final.{pair}.fq.gz",
            sample=SAMPLES_PE,
            pair="pe_pe pe_se".split()
        ),
        se_files = expand(
            QC + "{sample}.final.se.fq.gz",
            sample=SAMPLES_SE
        )


rule qc_fastqc_sample_pair:
    """
    Do FASTQC reports
    Uses --nogroup!
    One thread per fastq.gz file
    """
    input:
        fastq = QC + "{sample}.final.{pair}.fq.gz"
    output:
        zip = protected(QC + "{sample}.final.{pair}_fastqc.zip"),
        html = protected(QC + "{sample}.final.{pair}_fastqc.html")
    threads:
        1
    params:
        outdir = QC
    log:
        QC + "fastqc_{sample}_{pair}.log"
    benchmark:
        QC + "fastqc_{sample}_{pair}.bmk"
    conda:
        "qc.yml"
    shell:
        "fastqc --nogroup --outdir {params.outdir} {input.fastq} > {log} 2>&1"


rule qc_multiqc:
    input:
        pe_files = expand(
            QC + "{sample}.final.{pair}_fastqc.{extension}",
            sample=SAMPLES_PE,
            pair="pe_pe pe_se".split(),
            extension="html zip".split()
        ),
        se_files = expand(
            QC + "{sample}.final.se_fastqc.{extension}",
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
        "qc.yml"
    shell:
        "multiqc --title QC --filename {output.html} {params.folder} 2> {log}"


rule qc_doc:
    input:
        html = QC + "multiqc_report.html"


rule qc:
    """qc_results + qc_doc"""
    input:
        pe_files = expand(
            QC + "{sample}.final.{pair}_fastqc.{extension}",
            sample=SAMPLES_PE,
            pair="pe_pe pe_se".split(),
            extension="html zip".split()
        ),
        se_files = expand(
            QC + "{sample}.final.se.fq.gz",
            sample=SAMPLES_SE
        ),
        report = QC + "multiqc_report.html"
