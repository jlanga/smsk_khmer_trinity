rule diginorm_load_into_counting:
    """
    Build the hash table data structure from all the trimmed reads.
    Caution!: The --help says that it can be multithreaded but it raises
    errors!
    """
    input:
        fastqs = expand(
            QC + "{sample}.final.{pair}.fq.gz",
            sample=SAMPLES_PE,
            pair=PAIRS
        ) + expand(
            QC + "{sample}.final.se.fq.gz",
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


rule diginorm_normalize_by_median_sample_pe_pe:
    """
    Normalizes by median EACH FILE.
    Therefore one loads once per file the hash table.
    """
    input:
        fastq = QC + "{sample}.final.pe_pe.fq.gz",
        table = NORM + "diginorm_table.kh"
    output:
        fastq = temp(NORM + "{sample}.keep.pe_pe.fq.gz")
    threads:
        ALL_THREADS
    params:
        cutoff = config["diginorm_params"]["cutoff"],
        ksize = config["diginorm_params"]["ksize"],
        n_tables = config["diginorm_params"]["n_tables"],
        max_table_size = config["diginorm_params"]["max_table_size"],
    log:
        NORM + "normalize_by_median_{sample}.pe_pe.log"
    benchmark:
        NORM + "normalize_by_median_{sample}.pe_pe.bmk"
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


rule diginorm_normalize_by_median_sample_pe_se:
    """
    Normalizes by median EACH FILE.
    Therefore one loads once per file the hash table.
    """
    input:
        fastq = QC + "{sample}.final.pe_se.fq.gz",
        table = NORM + "diginorm_table.kh"
    output:
        fastq = temp(NORM + "{sample}.keep.pe_se.fq.gz")
    threads:
        ALL_THREADS
    params:
        cutoff = config["diginorm_params"]["cutoff"],
        ksize = config["diginorm_params"]["ksize"],
        n_tables = config["diginorm_params"]["n_tables"],
        max_table_size = config["diginorm_params"]["max_table_size"]
    log:
        NORM + "normalize_by_median_{sample}_pe_se.log"
    benchmark:
        NORM + "normalize_by_median_{sample}_pe_se.bmk"
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
        fastq = QC + "{sample}.final.se.fq.gz",
        table = NORM + "diginorm_table.kh"
    output:
        fastq = temp(NORM + "{sample}.keep.se.fq.gz")
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
        fastq = NORM + "{sample}.keep.{pair}.fq.gz",
        table = NORM + "diginorm_table.kh"
    output:
        fastq = temp(NORM + "{sample}.abundfilt.{pair}.fq.gz")
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
        fastq = NORM + "{sample}.abundfilt.pe_pe.fq.gz"
    output:
        fastq_pe = protected(NORM + "{sample}.final.pe_pe.fq.gz"),
        fastq_se = temp(NORM + "{sample}.temp.pe_se.fq.gz")
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
        from_norm = NORM + "{sample}.abundfilt.pe_se.fq.gz",
        from_paired = NORM + "{sample}.temp.pe_se.fq.gz"
    output:
        fastq = protected(NORM + "{sample}.final.pe_se.fq.gz")
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
        single = NORM + "{sample}.abundfilt.se.fq.gz"
    output:
        single = protected(NORM + "{sample}.final.se.fq.gz")
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
            NORM + "{sample}.final.{pair}.fq.gz",
            sample=SAMPLES_PE,
            pair=PAIRS
        ),
        expand(
            NORM + "{sample}.final.se.fq.gz",
            sample=SAMPLES_SE
        )


rule diginorm_fastqc_sample_pair:
    """
    Do FASTQC reports
    Uses --nogroup!
    One thread per fastq.gz file
    """
    input:
        fastq = NORM + "{sample}.final.{pair}.fq.gz"
    output:
        zip = protected(NORM + "{sample}.final.{pair}_fastqc.zip"),
        html = protected(NORM + "{sample}.final.{pair}_fastqc.html")
    log:
        NORM + "fastqc_{sample}_{pair}.log"
    benchmark:
        NORM + "fastqc_{sample}_{pair}.bmk"
    conda:
        "diginorm.yml"
    shell:
        "fastqc --nogroup --outdir {NORM} {input.fastq} > {log} 2>&1"


rule diginorm_multiqc:
    input:
        files_pe = expand(
            NORM + "{sample}.final.{pair}_fastqc.{extension}",
            sample=SAMPLES_PE,
            pair="pe_pe pe_se".split(),
            extension="html zip".split()
        ),
        files_se = expand(
            NORM + "{sample}.final.se_fastqc.{extension}",
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
        "diginorm.yml"
    shell:
        """
        multiqc \
            --title Diginorm \
            --filename {output.html} \
            {NORM} \
        2> {log}
        """


rule diginorm_doc:
    input:
        html = NORM + "multiqc_report.html"


rule diginorm:
    '''diginorm_results + diginorm_doc'''
    input:
        files_pe = expand(
            NORM + "{sample}.final.{pair}.fq.gz",
            sample=SAMPLES_PE,
            pair=PAIRS
        ),
        files_se = expand(
            NORM + "{sample}.final.se.fq.gz",
            sample=SAMPLES_SE
        ),
        html = NORM + "multiqc_report.html"
