rule diginorm_load_into_counting:
    """
    Build the hash table data structure from all the trimmed reads.
    Caution!: The --help says that it can be multithreaded but it raises
    errors!
    """
    input:
        fastqs = expand(
            QC + "{sample}_{pair}.fq.gz",
            sample=SAMPLES_PE,
            pair=PAIRS
        ) + expand(
            QC + "{sample}_se.fq.gz",
            sample=SAMPLES_SE
        )
    output:
        table = temp(NORM + "diginorm_table.kh"),
        info = temp(NORM + "diginorm_table.kh.info")
    threads:
        ALL_THREADS
    log:
        NORM + "load_into_counting.log"
    benchmark:
        NORM + "load_into_counting.bmk"
    params:
        ksize = config["diginorm_params"]["ksize"],
        max_table_size = config["diginorm_params"]["max_table_size"],
        n_tables = config["diginorm_params"]["n_tables"]
    conda:
        "diginorm.yml"
    shell:
        """
        load-into-counting.py \
            --ksize {params.ksize} \
            --n_tables {params.n_tables} \
            --max-tablesize {params.max_table_size} \
            --no-bigcount \
            --threads {threads} \
            {output.table} \
            {input.fastqs} \
        2> {log}
        """


rule diginorm_normalize_by_median_sample_pepe:
    """
    Normalizes by median EACH FILE.
    Therefore one loads once per file the hash table.
    """
    input:
        fastq = QC + "{sample}_pepe.fq.gz",
        table = NORM + "diginorm_table.kh"
    output:
        fastq = temp(NORM + "{sample}_pepe.keep.fq.gz")
    threads:
        ALL_THREADS
    params:
        cutoff = config["diginorm_params"]["cutoff"],
        ksize = config["diginorm_params"]["ksize"],
        n_tables = config["diginorm_params"]["n_tables"],
        max_table_size = config["diginorm_params"]["max_table_size"],
    log:
        NORM + "normalize_by_median_{sample}.pepe.log"
    benchmark:
        NORM + "normalize_by_median_{sample}.pepe.bmk"
    conda:
        "diginorm.yml"
    shell:
        """
        normalize-by-median.py \
            --ksize {params.ksize} \
            --n_tables {params.n_tables} \
            --max-tablesize {params.max_table_size} \
            --cutoff {params.cutoff} \
            --paired \
            --loadgraph {input.table} \
            --output >(pigz --best > {output.fastq}) \
            <(gzip --decompress --stdout {input.fastq}) \
        > {log} 2>&1
        """


rule diginorm_normalize_by_median_sample_pese:
    """
    Normalizes by median EACH FILE.
    Therefore one loads once per file the hash table.
    """
    input:
        fastq = QC + "{sample}_pese.fq.gz",
        table = NORM + "diginorm_table.kh"
    output:
        fastq = temp(NORM + "{sample}_pese.keep.fq.gz")
    threads:
        ALL_THREADS
    params:
        cutoff = config["diginorm_params"]["cutoff"],
        ksize = config["diginorm_params"]["ksize"],
        n_tables = config["diginorm_params"]["n_tables"],
        max_table_size = config["diginorm_params"]["max_table_size"]
    log:
        NORM + "normalize_by_median_{sample}_pese.log"
    benchmark:
        NORM + "normalize_by_median_{sample}_pese.bmk"
    conda:
        "diginorm.yml"
    shell:
        """
        normalize-by-median.py \
            --ksize {params.ksize} \
            --n_tables {params.n_tables} \
            --max-tablesize {params.max_table_size} \
            --cutoff {params.cutoff} \
            --loadgraph {input.table} \
            --output >(pigz --best > {output.fastq}) \
            <(gzip --decompress --stdout {input.fastq}) \
        > {log} 2>&1
        """


rule diginorm_normalize_by_median_sample_se:
    """
    Normalizes by median EACH FILE.
    Therefore one loads once per file the hash table.
    """
    input:
        fastq = QC + "{sample}_se.fq.gz",
        table = NORM + "diginorm_table.kh"
    output:
        fastq = temp(NORM + "{sample}_se.keep.fq.gz")
    threads:
        ALL_THREADS
    params:
        cutoff = config["diginorm_params"]["cutoff"],
        n_tables = config["diginorm_params"]["n_tables"],
        max_table_size = config["diginorm_params"]["max_table_size"],
        ksize = config["diginorm_params"]["ksize"]
    log:
        NORM + "normalize_by_median_{sample}_se.log"
    benchmark:
        NORM + "normalize_by_median_{sample}_se.bmk"
    conda:
        "diginorm.yml"
    shell:
        """
        normalize-by-median.py \
            --ksize {params.ksize} \
            --n_tables {params.n_tables} \
            --max-tablesize {params.max_table_size} \
            --cutoff {params.cutoff} \
            --loadgraph {input.table} \
            --output >(pigz --best > {output.fastq}) \
            <(gzip --decompress --stdout {input.fastq}) \
        > {log} 2>&1
        """


rule diginorm_filter_abund_sample_pair:
    """
    Removes erroneus k-mers.
    """
    input:
        fastq = NORM + "{sample}_{pair}.keep.fq.gz",
        table = NORM + "diginorm_table.kh"
    output:
        fastq = temp(NORM + "{sample}_{pair}.abundfilt.fq.gz")
    threads:
        ALL_THREADS
    priority:
        50
    log:
        NORM + "filter_abund_{sample}_{pair}.log"
    benchmark:
        NORM + "filter_abunt_{sample}_{pair}.bmk"
    conda:
        "diginorm.yml"
    shell:
        """
        filter-abund.py \
            --variable-coverage \
            --threads {threads} \
            --output >(pigz --best > {output.fastq}) \
            {input.table} \
            {input.fastq} \
        > {log} 2>&1
        """


rule diginorm_extract_paired_reads_sample:
    """
    Split the filtered reads into PE and SE.
    """
    input:
        fastq = NORM + "{sample}_pepe.abundfilt.fq.gz"
    output:
        fastq_pe = protected(NORM + "{sample}_pepe.fq.gz"),
        fastq_se = temp(NORM + "{sample}_pese.temp.fq.gz")
    threads:
        1
    log:
        NORM + "extract_paired_reads_{sample}.log"
    benchmark:
        NORM + "extract_paired_reads_{sample}.bmk"
    conda:
        "diginorm.yml"
    shell:
        """
        extract-paired-reads.py \
            --output-paired >(pigz --best > {output.fastq_pe}) \
            --output-single >(pigz --best > {output.fastq_se}) \
            {input.fastq} \
        > {log} 2>&1
        """


rule diginorm_merge_pe_single_reads_sample:
    """
    Put together the SE reads from the same sample
    """
    input:
        from_norm = NORM + "{sample}_pese.abundfilt.fq.gz",
        from_paired = NORM + "{sample}_pese.temp.fq.gz"
    output:
        fastq = protected(NORM + "{sample}_pese.fq.gz")
    log:
        NORM + "merge_single_reads_{sample}.log"
    benchmark:
        NORM + "merge_single_reads_{sample}.bmk"
    conda:
        "diginorm.yml"
    shell:
        """
        (pigz --decompress --stdout \
            {input.from_norm} \
            {input.from_paired} \
        | pigz --best \
        > {output.fastq}) \
        2> {log}
        """


rule dignorm_get_former_se_reads_sample:
    """
    Move the result of diginorm_extract_paired_reads for true SE reads
    to their final position.
    """
    input:
        single = NORM + "{sample}_se.abundfilt.fq.gz"
    output:
        single = protected(NORM + "{sample}_se.fq.gz")
    log:
        NORM + "get_former_se_reads_{sample}.log"
    benchmark:
        NORM + "get_former_se_reads_{sample}.bmk"
    conda:
        "diginorm.yml"
    shell:
        "mv {input.single} {output.single}"


rule diginorm_results:
    input:
        expand(
            NORM + "{sample}_{pair}.fq.gz",
            sample=SAMPLES_PE,
            pair=PAIRS
        ),
        expand(
            NORM + "{sample}_se.fq.gz",
            sample=SAMPLES_SE
        )


rule diginorm_doc:
    input:
        html = NORM + "multiqc_report.html"


rule diginorm:
    '''diginorm_results + diginorm_doc'''
    input:
        rules.diginorm_results.input,
        rules.diginorm_doc.input
