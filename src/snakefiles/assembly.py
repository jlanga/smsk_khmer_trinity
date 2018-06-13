rule assembly_split_pe_files:
    """
    Split pe_pe files into _1 and _2.
    """
    input:
        fastq_pe = NORM + "{sample}.final.pe_pe.fq.gz"
    output:
        left = ASSEMBLY + "{sample}_1.fq.gz",
        right = ASSEMBLY + "{sample}_2.fq.gz"
    params:
        left = "{sample}.final.pe_pe.fq.gz.1",
        right = "{sample}.final.pe_pe.fq.gz.2"
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


rule assembly_merge_right_and_left:
    """
    Generate the left.fq and right.fq
    left.fq = /1 reads from PE + all SE reads
    right.fq = /2 reads form PE
    Forced the use of gzip because pigz doesn't have --force option.
        --force is required because /dev/null isn't a file, doesn't even
        have an end, and isn't compressed.
    Test if SAMPLES_PE is empty because gzip/pigz may because otherwise
        it may be waiting something from stdout.
    """
    input:
        forward = expand(
            ASSEMBLY + "{sample}_1.fq.gz",
            sample=SAMPLES_PE
        ) if SAMPLES_PE else ["/dev/null"],
        reverse = expand(
            ASSEMBLY + "{sample}_2.fq.gz",
            sample=SAMPLES_PE
        ) if SAMPLES_PE else ["/dev/null"],
        single = expand(  # pe_se
            NORM + "{sample}.final.pe_se.fq.gz",
            sample=SAMPLES_PE
        ) + expand(  # se
            NORM + "{sample}.final.se.fq.gz",
            sample=SAMPLES_SE
        )
    output:
        left = ASSEMBLY + "left.fq",
        right = ASSEMBLY + "right.fq",
        single = ASSEMBLY + "single.fq"
    log:
        ASSEMBLY + "merge_right_and_left.log"
    benchmark:
        ASSEMBLY + "merge_right_and_left.bmk"
    conda:
        "assembly.yml"
    shell:
        """
        gzip --decompress --stdout {input.forward} > {output.left} 2> {log}
        gzip --decompress --stdout {input.reverse} > {output.right} 2>> {log}
        gzip --decompress --stdout {input.single} > {output.single} 2>> {log}
        """


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
        memory = config["trinity_params"]["memory"],
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


rule assembly_index_trinity:
    """
    Create a samtool index for the assembly
    """
    input:
        fasta = ASSEMBLY + "Trinity.fasta"
    output:
        fai = ASSEMBLY + "Trinity.fasta.fai"
    log:
        ASSEMBLY + "index.log"
    benchmark:
        ASSEMBLY + "index.bmk"
    conda:
        "assembly.yml"
    shell:
        "samtools faidx {input.fasta} 2> {log} 1>&2"
