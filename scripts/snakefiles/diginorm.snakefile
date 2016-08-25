rule diginorm_load_into_counting:
    """
    Build the hash table data structure from all the trimmed reads.
    Caution!: The --help says that it can be multithreaded but it raises
    errors!
    """
    input:
        fastqs = expand(
            QC_DIR + "{sample}.final.{pair}.fq.gz",
            sample = SAMPLES_PE,
            pair = PAIRS
        ) + expand(
            QC_DIR + "{sample}.final.se.fq.gz",
            sample = SAMPLES_SE
        )
    output:
        table = temp(NORM_DIR + "diginorm_table.kh"),
        info  = temp(NORM_DIR + "diginorm_table.kh.info")
    threads:
        1
    log:
        NORM_DOC + "load_into_counting.log"
    benchmark:
        NORM_DOC + "load_into_counting.json"
    params:
        ksize=    config["diginorm_params"]["ksize"],
        hashsize= config["diginorm_params"]["hashsize"],
        n_hashes= config["diginorm_params"]["n_hashes"]
    shell:
        """
        load-into-counting.py \
            --ksize {params.ksize} \
            --n_hashes {params.n_hashes} \
            --hashsize {params.hashsize} \
            --no-bigcount \
            {output.table} \
            {input.fastqs} \
        >  {log} 2>&1
        """



rule diginorm_normalize_by_median_sample_pe_pe:
    """
    Normalizes by median EACH FILE.
    Therefore one loads once per file the hash table.
    """
    input:
        fastq = QC_DIR + "{sample}.final.pe_pe.fq.gz",
        table = NORM_DIR + "diginorm_table.kh"
    output:
        fastq = temp(NORM_DIR + "{sample}.keep.pe_pe.fq.gz")
    threads:
        BLOCK_THREADS # Block
    params:
        cutoff   = config["diginorm_params"]["cutoff"],
        ksize    = config["diginorm_params"]["ksize"],
        n_hashes = config["diginorm_params"]["n_hashes"],
        hashsize = config["diginorm_params"]["hashsize"],
        keep_fq  = "{sample}.final.pe_pe.fq.gz.keep"
    log:
        NORM_DOC + "normalize_by_median_{sample}.pe_pe.log"
    benchmark:
        NORM_DOC + "normalize_by_median_{sample}.pe_pe.json"
    shell:
        """
        normalize-by-median.py \
            --ksize {params.ksize} \
            --n_hashes {params.n_hashes} \
            --hashsize {params.hashsize} \
            --cutoff {params.cutoff} \
            --paired \
            --loadhash {input.table} \
            {input.fastq} \
        > {log} 2>&1

        pigz -9c \
            {params.keep_fq} \
        > {output.fastq} \
        2>> {log}

        rm {params.keep_fq}
        """



rule diginorm_normalize_by_median_sample_pe_se:
    """
    Normalizes by median EACH FILE.
    Therefore one loads once per file the hash table.
    """
    input:
        fastq = QC_DIR + "{sample}.final.pe_se.fq.gz",
        table = NORM_DIR + "diginorm_table.kh"
    output:
        fastq = temp(NORM_DIR + "{sample}.keep.pe_se.fq.gz")
    threads:
        BLOCK_THREADS # Block excessive RAM usage
    params:
        cutoff   = config["diginorm_params"]["cutoff"],
        ksize    = config["diginorm_params"]["ksize"],
        n_hashes = config["diginorm_params"]["n_hashes"],
        hashsize = config["diginorm_params"]["hashsize"],
        keep_fq  = "{sample}.final.pe_se.fq.gz.keep"
    log:
        NORM_DOC + "normalize_by_median_{sample}_pe_se.log"
    benchmark:
        NORM_DOC + "normalize_by_median_{sample}_pe_se.json"
    shell:
        """
        normalize-by-median.py \
            --ksize {params.ksize} \
            --n_hashes {params.n_hashes} \
            --hashsize {params.hashsize} \
            --cutoff {params.cutoff} \
            --loadhash {input.table} \
            {input.fastq} \
        > {log} 2>&1

        pigz -9c \
            {params.keep_fq} \
        > {output.fastq} \
        2>> {log}

        rm {params.keep_fq}
        """



rule diginorm_normalize_by_median_sample_se:
    """
    Normalizes by median EACH FILE.
    Therefore one loads once per file the hash table.
    """
    input:
        fastq = QC_DIR + "{sample}.final.se.fq.gz",
        table = NORM_DIR + "diginorm_table.kh"
    output:
        fastq = temp(NORM_DIR + "{sample}.keep.se.fq.gz")
    threads:
        BLOCK_THREADS # Block excessive RAM usage
    params:
        cutoff   = config["diginorm_params"]["cutoff"],
        ksize    = config["diginorm_params"]["ksize"],
        n_hashes = config["diginorm_params"]["n_hashes"],
        hashsize = config["diginorm_params"]["hashsize"],
        keep_fq  = "{sample}.final.se.fq.gz.keep"
    log:
        NORM_DOC + "normalize_by_median_{sample}_se.log"
    benchmark:
        NORM_DOC + "normalize_by_median_{sample}_se.json"
    shell:
        """
        normalize-by-median.py \
            --ksize {params.ksize} \
            --n_hashes {params.n_hashes} \
            --hashsize {params.hashsize} \
            --cutoff {params.cutoff} \
            --loadhash {input.table} \
            {input.fastq} \
        > {log} 2>&1

        pigz -9c \
            {params.keep_fq} \
        > {output.fastq} \
        2>> {log}

        rm {params.keep_fq}
        """



rule diginorm_filter_abund_sample_pair:
    """
    Removes erroneus k-mers.
    """
    input:
        fastq = NORM_DIR + "{sample}.keep.{pair}.fq.gz",
        table = NORM_DIR + "diginorm_table.kh"
    output:
        fastq = temp(NORM_DIR + "{sample}.abundfilt.{pair}.fq.gz")
    threads:
        BLOCK_THREADS # BLOCK
    params:
        abundfilt_fq = "{sample}.keep.{pair}.fq.gz.abundfilt"
    log:
        NORM_DOC + "filter_abund_{sample}_{pair}.log"
    benchmark:
        NORM_DOC + "filter_abunt_{sample}_{pair}.json"
    shell:
        """
        filter-abund.py \
            --variable-coverage \
            {input.table} \
            {input.fastq} \
        > {log} 2>&1

        pigz -9c \
            {params.abundfilt_fq} \
        > {output.fastq}  \
        2>> {log}

        rm {params.abundfilt_fq}
        """



rule diginorm_extract_paired_reads_sample:
    """
    Split the filtered reads into PE and SE.
    """
    input:
        fastq = NORM_DIR + "{sample}.abundfilt.pe_pe.fq.gz"
    output:
        fastq_pe = protected(NORM_DIR + "{sample}.final.pe_pe.fq.gz"),
        fastq_se = temp(NORM_DIR + "{sample}.temp.pe_se.fq.gz")
    threads:
        1
    params:
        fastq_pe = "{sample}.abundfilt.pe_pe.fq.gz.pe",
        fastq_se = "{sample}.abundfilt.pe_pe.fq.gz.se"
    log:
        NORM_DOC + "extract_paired_reads_{sample}.log"
    benchmark:
        NORM_DOC + "extract_paired_reads_{sample}.json"
    shell:
        """
        extract-paired-reads.py \
            {input.fastq} \
        > {log} 2>&1

        pigz -9c {params.fastq_pe} > {output.fastq_pe}
        pigz -9c {params.fastq_se} > {output.fastq_se}

        rm {params.fastq_pe} {params.fastq_se}
        """



rule diginorm_merge_pe_single_reads_sample:
    """
    Put together the SE reads from the same sample
    """
    input:
        from_norm=   NORM_DIR + "{sample}.abundfilt.pe_se.fq.gz",
        from_paired= NORM_DIR + "{sample}.temp.pe_se.fq.gz"
    output:
        fastq = protected(NORM_DIR + "{sample}.final.pe_se.fq.gz")
    threads:
        1
    log:
        NORM_DOC + "merge_single_reads_{sample}.log"
    benchmark:
        NORM_DOC + "merge_single_reads_{sample}.json"
    shell:
        """
        cp {input.from_norm} {output.fastq}
        pigz -dc {input.from_paired} |
        pigz -9 >> {output.fastq}
        """



rule dignorm_get_former_se_reads_sample:
    """
    Move the result of diginorm_extract_paired_reads for true SE reads
    to their final position.
    """
    input:
        single= NORM_DIR + "{sample}.abundfilt.se.fq.gz"
    output:
        single= protected(NORM_DIR + "{sample}.final.se.fq.gz")
    threads:
        1
    log:
        NORM_DOC + "get_former_se_reads_{sample}.log"
    benchmark:
        NORM_DOC + "get_former_se_reads_{sample}.json"
    shell:
        """
        mv {input.single} {output.single}
        """


rule diginorm_results:
    input:
        expand(
            NORM_DIR + "{sample}.final.{pair}.fq.gz",
            sample = SAMPLES_PE,
            pair = PAIRS
        ),
        expand(
            NORM_DIR + "{sample}.final.se.fq.gz",
            sample = SAMPLES_SE
        )



rule diginorm_fastqc_sample_pair:
    """
    Do FASTQC reports
    Uses --nogroup!
    One thread per fastq.gz file
    """
    input:
        fastq = NORM_DIR + "{sample}.final.{pair}.fq.gz"
    output:
        zip   = protected(NORM_DOC + "{sample}.final.{pair}_fastqc.zip"),
        html  = protected(NORM_DOC + "{sample}.final.{pair}_fastqc.html")
    threads:
        1
    params:
        outdir = NORM_DOC
    log:
        NORM_DOC + "fastqc_{sample}_{pair}.log"
    benchmark:
        NORM_DOC + "fastqc_{sample}_{pair}.json"
    shell:
        "fastqc "
            "--nogroup "
            "--outdir {params.outdir} "
            "{input.fastq} "
        "> {log} 2>&1"



rule diginorm_report:
    input:
        expand(
            NORM_DOC + "{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_PE,
            pair = "pe_pe pe_se".split(),
            extension = "html zip".split()
        ),
        expand(
            NORM_DOC + "{sample}.final.se_fastqc.{extension}",
            sample = SAMPLES_SE,
            extension = "html zip".split()
        )



rule diginorm:
    input:
        expand(
            NORM_DIR + "{sample}.final.{pair}.fq.gz",
            sample = SAMPLES_PE,
            pair = PAIRS
        ),
        expand(
            NORM_DIR + "{sample}.final.se.fq.gz",
            sample = SAMPLES_SE
        ),
        expand(
            NORM_DOC + "{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_PE,
            pair = "pe_pe pe_se".split(),
            extension = "html zip".split()
        ),
        expand(
            NORM_DOC + "{sample}.final.se_fastqc.{extension}",
            sample = SAMPLES_SE,
            extension = "html zip".split()
        )
