rule assembly_split_pe_files:
    """
    Split pepe files into _1 and _2.
    """
    input:
        fastq_pe = NORM + "{sample}_pepe.fq.gz"
    output:
        left = ASSEMBLY + "{sample}_1.fq.gz",
        right = ASSEMBLY + "{sample}_2.fq.gz"
    log:
        ASSEMBLY + "split_pe_files_{sample}.log"
    benchmark:
        ASSEMBLY + "split_pe_files_{sample}.bmk"
    conda:
        "assembly.yml"
    shell:
        """
        split-paired-reads.py \
            --output-first >(pigz --best > {output.left}) \
            --output-second >(pigz --best > {output.right}) \
            {input.fastq_pe} \
        > {log} 2>&1
        sleep 5
        """


rule assembly_merge_left:
    input:
        expand(
            ASSEMBLY + "{sample}_1.fq.gz",
            sample=SAMPLES_PE
        ) if SAMPLES_PE else ["/dev/null"]
    output: ASSEMBLY + "left.fq"
    log: ASSEMBLY + "merge_left.log"
    benchmark: ASSEMBLY + "merge_left.bmk"
    shell: "gzip --decompress --stdout {input} > {output} 2> {log}"


rule assembly_merge_right:
    input:
        expand(
            ASSEMBLY + "{sample}_2.fq.gz",
            sample=SAMPLES_PE
        ) if SAMPLES_PE else ["/dev/null"]
    output: ASSEMBLY + "right.fq"
    log: ASSEMBLY + "merge_right.log"
    benchmark: ASSEMBLY + "merge_right.bmk"
    shell: "gzip --decompress --stdout {input} > {output} 2> {log}"


rule assembly_merge_single:
    input:
        expand(
            ASSEMBLY + "{sample}_2.fq.gz",
            sample=SAMPLES_PE
        ) if SAMPLES_PE else ["/dev/null"]
    output: ASSEMBLY + "single.fq"
    log: ASSEMBLY + "merge_single.log"
    benchmark: ASSEMBLY + "merge_single.bmk"
    shell: "gzip --decompress --stdout {input} > {output} 2> {log}"


rule assembly_run_trinity:
    """
    Assembly reads with Trinity.
    Notes on hardcoded settings:
        - Runs on paired end mode
        - Expect fastq files as inputs (left and right)
        - Does the full cleanup so it only remains a fasta file.
    """
    input:
        left = ASSEMBLY + "left.fq",
        right = ASSEMBLY + "right.fq",
        # single = ASSEMBLY + "single.fq"
    output:
        fasta = protected(ASSEMBLY + "Trinity.fasta")
    threads:
        ALL_THREADS
    priority:
        50
    params:
        memory = params["trinity"]["memory"],
        outdir = ASSEMBLY + "trinity_out_dir"
    log:
        ASSEMBLY + "run_trinity.log"
    benchmark:
        ASSEMBLY + "run_trinity.bmk"
    conda:
        "assembly.yml"
    shell:
        """
        Trinity \
            --seqType fq \
            --no_normalize_reads \
            --max_memory {params.memory} \
            --left {input.left} \
            --right {input.right} \
            --CPU {threads} \
            --full_cleanup \
            --output {params.outdir} \
        > {log}
        mv {params.outdir}.Trinity.fasta {output.fasta}
        """


rule assembly_gene_to_trans_map:
    """
    Create the gene_id TAB transcript_id file
    """
    input:
        fasta = ASSEMBLY + "Trinity.fasta"
    output:
        tsv = ASSEMBLY + "Trinity_gene_to_trans.tsv"
    log:
        ASSEMBLY + "gene_to_trans_map.log"
    benchmark:
        ASSEMBLY + "gene_to_trans_map.bmk"
    conda:
        "assembly.yml"
    shell:
        "get_Trinity_gene_to_trans_map.pl < {input.fasta} > {output.tsv} "
        "2> {log}"
