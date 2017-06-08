rule diginorm_load_into_counting:
    """
    Build the hash table data structure from all the trimmed reads.
    Caution!: The --help says that it can be multithreaded but it raises
    errors!
    """
    input:
        fastqs = expand(
            qc + "{sample}.final.{pair}.fq.gz",
            sample = SAMPLES_PE,
            pair = PAIRS
        ) + expand(
            qc + "{sample}.final.se.fq.gz",
            sample = SAMPLES_SE
        )
    output:
        table = temp(norm + "diginorm_table.kh"),
        info  = temp(norm + "diginorm_table.kh.info")
    threads:
        ALL_THREADS
    priority:
        50
    log:
        norm + "load_into_counting.log"
    benchmark:
        norm + "load_into_counting.json"
    params:
        ksize= config["diginorm_params"]["ksize"],
        max_table_size= config["diginorm_params"]["max_table_size"],
        n_tables= config["diginorm_params"]["n_tables"]
    shell:
        "load-into-counting.py "
            "--ksize {params.ksize} "
            "--n_tables {params.n_tables} "
            "--max-tablesize {params.max_table_size} "
            "--no-bigcount "
            "--threads {threads} "
            "{output.table} "
            "{input.fastqs} "
        "> {log} 2>&1"



rule diginorm_normalize_by_median_sample_pe_pe:
    """
    Normalizes by median EACH FILE.
    Therefore one loads once per file the hash table.
    """
    input:
        fastq = qc + "{sample}.final.pe_pe.fq.gz",
        table = norm + "diginorm_table.kh"
    output:
        fastq = temp(norm + "{sample}.keep.pe_pe.fq.gz")
    threads:
        ALL_THREADS
    params:
        cutoff   = config["diginorm_params"]["cutoff"],
        ksize    = config["diginorm_params"]["ksize"],
        n_tables = config["diginorm_params"]["n_tables"],
        max_table_size = config["diginorm_params"]["max_table_size"],
    log:
        norm + "normalize_by_median_{sample}.pe_pe.log"
    benchmark:
        norm + "normalize_by_median_{sample}.pe_pe.json"
    shell:
        "normalize-by-median.py "
            "--ksize {params.ksize} "
            "--n_tables {params.n_tables} "
            "--max-tablesize {params.max_table_size} "
            "--cutoff {params.cutoff} "
            "--paired "
            "--loadgraph {input.table} "
            "--output >(pigz --best > {output.fastq}) "
            "<(gzip --decompress --stdout {input.fastq}) "
        "> {log} 2>&1 "



rule diginorm_normalize_by_median_sample_pe_se:
    """
    Normalizes by median EACH FILE.
    Therefore one loads once per file the hash table.
    """
    input:
        fastq = qc + "{sample}.final.pe_se.fq.gz",
        table = norm + "diginorm_table.kh"
    output:
        fastq = temp(norm + "{sample}.keep.pe_se.fq.gz")
    threads:
        ALL_THREADS
    params:
        cutoff   = config["diginorm_params"]["cutoff"],
        ksize    = config["diginorm_params"]["ksize"],
        n_tables = config["diginorm_params"]["n_tables"],
        max_table_size = config["diginorm_params"]["max_table_size"]
    log:
        norm + "normalize_by_median_{sample}_pe_se.log"
    benchmark:
        norm + "normalize_by_median_{sample}_pe_se.json"
    shell:
        "normalize-by-median.py "
            "--ksize {params.ksize} "
            "--n_tables {params.n_tables} "
            "--max-tablesize {params.max_table_size} "
            "--cutoff {params.cutoff} "
            "--loadgraph {input.table} "
            "--output >(pigz --best > {output.fastq}) "
            "<(gzip --decompress --stdout {input.fastq}) "
        "> {log} 2>&1"



rule diginorm_normalize_by_median_sample_se:
    """
    Normalizes by median EACH FILE.
    Therefore one loads once per file the hash table.
    """
    input:
        fastq = qc + "{sample}.final.se.fq.gz",
        table = norm + "diginorm_table.kh"
    output:
        fastq = temp(norm + "{sample}.keep.se.fq.gz")
    threads:
        ALL_THREADS
    params:
        cutoff   = config["diginorm_params"]["cutoff"],
        ksize    = config["diginorm_params"]["ksize"],
        n_tables = config["diginorm_params"]["n_tables"],
        max_table_size = config["diginorm_params"]["max_table_size"]
    log:
        norm + "normalize_by_median_{sample}_se.log"
    benchmark:
        norm + "normalize_by_median_{sample}_se.json"
    shell:
        "normalize-by-median.py "
            "--ksize {params.ksize} "
            "--n_tables {params.n_tables} "
            "--max-tablesize {params.max_table_size} "
            "--cutoff {params.cutoff} "
            "--loadgraph {input.table} "
            "--output >(pigz --best > {output.fastq}) "
            "<(gzip --decompress --stdout {input.fastq}) "
        "> {log} 2>&1"




rule diginorm_filter_abund_sample_pair:
    """
    Removes erroneus k-mers.
    """
    input:
        fastq = norm + "{sample}.keep.{pair}.fq.gz",
        table = norm + "diginorm_table.kh"
    output:
        fastq = temp(norm + "{sample}.abundfilt.{pair}.fq.gz")
    threads:
        ALL_THREADS
    priority:
        50
    log:
        norm + "filter_abund_{sample}_{pair}.log"
    benchmark:
        norm + "filter_abunt_{sample}_{pair}.json"
    shell:
        "filter-abund.py "
            "--variable-coverage "
            "--threads {threads} "
            "--output >(pigz --best > {output.fastq}) "
            "{input.table} "
            "{input.fastq} "
        "> {log} 2>&1"



rule diginorm_extract_paired_reads_sample:
    """
    Split the filtered reads into PE and SE.
    """
    input:
        fastq = norm + "{sample}.abundfilt.pe_pe.fq.gz"
    output:
        fastq_pe = protected(norm + "{sample}.final.pe_pe.fq.gz"),
        fastq_se = temp(norm + "{sample}.temp.pe_se.fq.gz")
    threads:
        2
    log:
        norm + "extract_paired_reads_{sample}.log"
    benchmark:
        norm + "extract_paired_reads_{sample}.json"
    shell:
        "extract-paired-reads.py "
            "--output-paired >(pigz --best > {output.fastq_pe}) "
            "--output-single >(pigz --best > {output.fastq_se}) "
            "{input.fastq} "
        "> {log} 2>&1"


rule diginorm_merge_pe_single_reads_sample:
    """
    Put together the SE reads from the same sample
    """
    input:
        from_norm=   norm + "{sample}.abundfilt.pe_se.fq.gz",
        from_paired= norm + "{sample}.temp.pe_se.fq.gz"
    output:
        fastq = protected(norm + "{sample}.final.pe_se.fq.gz")
    log:
        norm + "merge_single_reads_{sample}.log"
    benchmark:
        norm + "merge_single_reads_{sample}.json"
    shell:
        "(pigz --decompress --stdout "
            "{input.from_norm} "
            "{input.from_paired} "
        "| pigz --best "
        "> {output.fastq}) "
        "2> {log}"



rule dignorm_get_former_se_reads_sample:
    """
    Move the result of diginorm_extract_paired_reads for true SE reads
    to their final position.
    """
    input:
        single= norm + "{sample}.abundfilt.se.fq.gz"
    output:
        single= protected(norm + "{sample}.final.se.fq.gz")
    log:
        norm + "get_former_se_reads_{sample}.log"
    benchmark:
        norm + "get_former_se_reads_{sample}.json"
    shell:
        "mv {input.single} {output.single}"



rule diginorm_results:
    input:
        expand(
            norm + "{sample}.final.{pair}.fq.gz",
            sample = SAMPLES_PE,
            pair = PAIRS
        ),
        expand(
            norm + "{sample}.final.se.fq.gz",
            sample = SAMPLES_SE
        )



rule diginorm_fastqc_sample_pair:
    """
    Do FASTQC reports
    Uses --nogroup!
    One thread per fastq.gz file
    """
    input:
        fastq = norm + "{sample}.final.{pair}.fq.gz"
    output:
        zip   = protected(norm + "{sample}.final.{pair}_fastqc.zip"),
        html  = protected(norm + "{sample}.final.{pair}_fastqc.html")
    params:
        outdir = norm
    log:
        norm + "fastqc_{sample}_{pair}.log"
    benchmark:
        norm + "fastqc_{sample}_{pair}.json"
    shell:
        "fastqc "
            "--nogroup "
            "--outdir {params.outdir} "
            "{input.fastq} "
        "> {log} 2>&1"



rule diginorm_multiqc:
    input:
        files_pe = expand(
            norm + "{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_PE,
            pair = "pe_pe pe_se".split(),
            extension = "html zip".split()
        ),
        files_se = expand(
            norm + "{sample}.final.se_fastqc.{extension}",
            sample = SAMPLES_SE,
            extension = "html zip".split()
        )
    output:
        html= protected(norm + "multiqc_report.html")
    params:
        folder = norm
    log:
        norm + "multiqc.log"
    benchmark:
        norm + "multiqc.json"
    shell:
        "multiqc "
            "--title Diginorm "
            "--filename {output.html} "
            "{params.folder} "
        "2> {log}"



rule diginorm_doc:
    input:
        html= norm + "multiqc_report.html"



rule diginorm:
    '''diginorm_results + diginorm_doc'''
    input:
        files_pe = expand(
            norm + "{sample}.final.{pair}.fq.gz",
            sample = SAMPLES_PE,
            pair = PAIRS
        ),
        files_se = expand(
            norm + "{sample}.final.se.fq.gz",
            sample = SAMPLES_SE
        ),
        html= norm + "multiqc_report.html"
