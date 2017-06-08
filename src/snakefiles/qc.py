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
        forward = raw + "{sample}_1.fq.gz",
        reverse = raw + "{sample}_2.fq.gz"
    output:
        forward     = temp(qc + "{sample}_1.fq.gz"),
        reverse     = temp(qc + "{sample}_2.fq.gz"),
        unpaired    = protected(qc + "{sample}.final.pe_se.fq.gz")
    params:
        unpaired_1  = qc + "{sample}_3.fq.gz",
        unpaired_2  = qc + "{sample}_4.fq.gz",
        adaptor     = lambda wildcards: config["samples_pe"][wildcards.sample]["adaptor"],
        phred       = lambda wildcards: config["samples_pe"][wildcards.sample]["phred"],
        trimmomatic_params = config["trimmomatic_params"]
    log:
        qc + "trimmomatic_pe_{sample}.log"
    benchmark:
        qc + "trimmomatic_pe_{sample}.json"
    priority:
        50
    threads:
        ALL_THREADS
    shell:
        "(trimmomatic PE "
            "-threads {threads} "
            "-{params.phred} "
            "<(gzip --decompress --stdout {input.forward}) "
            "<(gzip --decompress --stdout {input.reverse}) "
            ">(cut --fields 1 --delimiter \" \" | pigz --best > {output.forward}) "
            "{params.unpaired_1} "
            ">(cut --fields 1 --delimiter \" \" | pigz --best > {output.reverse}) "
            "{params.unpaired_2} "
            "ILLUMINACLIP:{params.adaptor}:2:30:10 "
            "{params.trimmomatic_params} ; "
        "gzip --decompress --stdout {params.unpaired_1} {params.unpaired_2} "
            "| cut --fields 1 --delimiter \" \" "
            "| pigz --best "
        "> {output.unpaired} ; "
        "rm {params.unpaired_1} {params.unpaired_2}) "
        "2> {log}"



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
        single = raw + "{sample}_se.fq.gz",
    output:
        single = protected(qc + "{sample}.final.se.fq.gz")
    params:
        adaptor = lambda wildcards: config["samples_se"][wildcards.sample]["adaptor"],
        phred = lambda wildcards: config["samples_se"][wildcards.sample]["phred"],
        trimmomatic_params = config["trimmomatic_params"]
    priority:
        50
    log:
        qc + "trimmomatic_se_{sample}.log"
    benchmark:
        qc + "trimmomatic_se_{sample}.json"
    threads:
        ALL_THREADS
    shell:
        "trimmomatic SE "
            "-threads {threads} "
            "-{params.phred} "
            "<(gzip --decompress --stdout {input.single}) "
            ">(cut --fields 1 --delimiter \" \" | pigz --best > {output.single}) "
            "ILLUMINACLIP:{params.adaptor}:2:30:10 "
            "{params.trimmomatic_params} "
        "2> {log}"



rule qc_decompress_pe_sample:
    input:
        forward = qc + "{sample}_1.fq.gz",
        reverse = qc + "{sample}_2.fq.gz"
    output:
        forward = qc + "{sample}_1.fq",
        reverse = qc + "{sample}_2.fq"
    log:
        qc + "decompress_pe_{sample}.log"
    benchmark:
        qc + "decompress_pe_{sample}.json"
    shell:
        "gzip --decompress --keep {input.forward} {input.reverse} 2> {log}"



rule qc_interleave_pe_pe:
    """
    From the adaptor free _1 and _2 , interleave the reads.
    Read the inputs, interleave, filter the stream and compress.
    """
    input:
        forward= qc + "{sample}_1.fq.gz",
        reverse= qc + "{sample}_2.fq.gz"
    output:
        interleaved= protected(qc + "{sample}.final.pe_pe.fq.gz")
    threads:
        4
    log:
        qc + "interleave_pe_{sample}.log"
    benchmark:
        qc + "interleave_pe_{sample}.json"
    shell:
        "(interleave-reads.py "
            "<(gzip --decompress --stdout {input.forward}) "
            "<(gzip --decompress --stdout {input.reverse}) "
            "| pigz --best "
        "> {output.interleaved}) "
        "2> {log}"



rule qc_results:
    '''Generate only the resulting files, not the reports'''
    input:
        pe_files = expand(
            qc + "{sample}.final.{pair}.fq.gz",
            sample = SAMPLES_PE,
            pair = "pe_pe pe_se".split()
        ),
        se_files = expand(
            qc + "{sample}.final.se.fq.gz",
            sample = SAMPLES_SE
        )



rule qc_fastqc_sample_pair:
    """
    Do FASTQC reports
    Uses --nogroup!
    One thread per fastq.gz file
    """
    input:
        fastq = qc + "{sample}.final.{pair}.fq.gz"
    output:
        zip   = protected(qc + "{sample}.final.{pair}_fastqc.zip"),
        html  = protected(qc + "{sample}.final.{pair}_fastqc.html")
    threads:
        1
    params:
        outdir = qc
    log:
        qc + "fastqc_{sample}_{pair}.log"
    benchmark:
        qc + "fastqc_{sample}_{pair}.json"
    shell:
        "fastqc "
            "--nogroup "
            "--outdir {params.outdir} "
            "{input.fastq} "
        "> {log} 2>&1"



rule qc_multiqc:
    input:
        pe_files = expand(
            qc + "{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_PE,
            pair = "pe_pe pe_se".split(),
            extension = "html zip".split()
        ),
        se_files = expand(
            qc + "{sample}.final.se_fastqc.{extension}",
            sample = SAMPLES_SE,
            extension = "html zip".split()
        )
    output:
        html= protected(qc + "multiqc_report.html")
    threads:
        1
    params:
        folder = qc
    log:
        qc + "multiqc.log"
    benchmark:
        qc + "multiqc.json"
    shell:
        "multiqc "
            "--title QC "
            "--filename {output.html} "
            "{params.folder} "
        "2> {log}"



rule qc_doc:
    input:
        html= qc + "multiqc_report.html"



rule qc:
    """qc_results + qc_doc"""
    input:
        pe_files = expand(
            qc + "{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_PE,
            pair = "pe_pe pe_se".split(),
            extension = "html zip".split()
        ),
        se_files = expand(
            qc + "{sample}.final.se.fq.gz",
            sample = SAMPLES_SE
        ),
        report = qc + "multiqc_report.html"
