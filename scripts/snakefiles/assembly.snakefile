rule assembly_split_pe_files:
    """
    Split pe_pe files into _1 and _2.
    """
    input:
        fastq_pe = NORM_DIR + "{sample}.final.pe_pe.fq.gz"
    output:
        left  = temp(ASSEMBLY_DIR + "{sample}_1.fq.gz"),
        right = temp(ASSEMBLY_DIR + "{sample}_2.fq.gz")
    threads:
        1
    params:
        left  = "{sample}.final.pe_pe.fq.gz.1",
        right = "{sample}.final.pe_pe.fq.gz.2"
    log:
        ASSEMBLY_DOC + "split_pe_files_{sample}.log"
    benchmark:
        ASSEMBLY_DOC + "split_pe_files_{sample}.json"
    shell:
        """
        split-paired-reads.py \
            {input.fastq_pe} \
        > {log} 2>&1

        pigz -9c {params.left}  > {output.left}
        pigz -9c {params.right} > {output.right}

        rm {params.left} {params.right}
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
        forward=  expand(
            ASSEMBLY_DIR + "{sample}_1.fq.gz",
            sample = SAMPLES_PE
        ) if SAMPLES_PE  else ["/dev/null"],
        reverse=  expand(
            ASSEMBLY_DIR + "{sample}_2.fq.gz",
            sample = SAMPLES_PE
        ) if SAMPLES_PE  else ["/dev/null"],
        unpaired= expand( # pe_se
            NORM_DIR + "{sample}.final.pe_se.fq.gz",
            sample = SAMPLES_PE
        ) + expand( # se
            NORM_DIR + "{sample}.final.se.fq.gz",
            sample = SAMPLES_SE
        )
    output:
        left=  temp(ASSEMBLY_DIR + "left.fq"),
        right= temp(ASSEMBLY_DIR + "right.fq")
    threads:
        1
    log:
        ASSEMBLY_DOC + "merge_right_and_left.log"
    benchmark:
        ASSEMBLY_DOC + "merge_right_and_left.json"
    shell:
        """
        gzip -dfc {input.forward} {input.unpaired} > {output.left} 2> {log}
        gzip -dfc {input.reverse} > {output.right} 2>> {log}
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
        left  = ASSEMBLY_DIR + "left.fq",
        right = ASSEMBLY_DIR + "right.fq"
    output:
        fasta = protected(ASSEMBLY_DIR + "Trinity.fasta")
    threads:
        ALL_THREADS
    params:
        memory= config["trinity_params"]["memory"],
        outdir= ASSEMBLY_DIR + "trinity_out_dir"
    log:
        ASSEMBLY_DOC + "run_trinity.log"
    benchmark:
        ASSEMBLY_DOC + "run_trinity.json"
    shell:
        """
        ./bin/Trinity \
            --seqType fq \
            --max_memory {params.memory} \
            --left {input.left} \
            --right {input.right} \
            --CPU {threads} \
            --full_cleanup \
            --output {params.outdir} \
        > {log}

        mv {params.outdir}.Trinity.fasta {output.fasta}

        rm {input.left}.readcount {input.right}.readcount
        """
